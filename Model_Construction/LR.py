import pandas as pd
from sklearn.model_selection import train_test_split, RandomizedSearchCV, StratifiedKFold
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import SelectKBest, f_classif, RFE, SelectFromModel
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, classification_report, roc_curve, auc
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import randint, uniform
import joblib  # 用于保存和加载模型
import os
from lifelines import KaplanMeierFitter, CoxPHFitter

# 设置中文字体
plt.rcParams['font.sans-serif'] = ['SimHei']
plt.rcParams['axes.unicode_minus'] = False

# 1. 读取训练数据
df_train = pd.read_excel("test.xlsx")

# 2. 数据预处理
def preprocess_data(df):
    df['Group'] = df['Group'].apply(lambda x: 1 if x == 'Cluster3' else 0)
    X = df.drop(['ID', 'Group'], axis=1)
    y = df['Group']
    return X, y

X, y = preprocess_data(df_train)

# 3. 划分训练集和验证集 (7:3) - 用于参数优化和评估
X_train, X_val, y_train, y_val = train_test_split(X, y, test_size=0.3, random_state=42)

# 4. 记录训练数据的列名
train_cols = X_train.columns

# 5. 定义特征选择方法
feature_selection_methods = {
    'SelectKBest': SelectKBest(f_classif, k=20),  # 选择 20 个最佳特征
    'RFE_LR': RFE(estimator=LogisticRegression(max_iter=1000), n_features_to_select=20),  # 使用逻辑回归进行递归特征消除
    'SelectFromModel_RF': SelectFromModel(RandomForestClassifier(n_estimators=500, random_state=42)),  # 基于随机森林模型的重要性进行特征选择
    'SelectFromModel_LR': SelectFromModel(LogisticRegression(penalty='l1', solver='liblinear', max_iter=1000)),  # 基于逻辑回归模型的 L1 正则化进行特征选择
    'RFE_RF': RFE(estimator=RandomForestClassifier(n_estimators=500, random_state=42), n_features_to_select=20),  # 使用随机森林进行递归特征消除
}

# 6. 训练和评估模型，并预测外部数据
for method_name, method in feature_selection_methods.items():
    print(f"使用特征选择方法: {method_name}")

    # 对训练数据进行特征选择
    X_train_selected = method.fit_transform(X_train, y_train)
    X_val_selected = method.transform(X_val)
    selected_feature_names = X_train.columns[method.get_support()]
    print(f"选择的特征: {selected_feature_names}")

    # 随机森林模型参数优化和交叉验证
    param_distributions = {
        'n_estimators': randint(100, 1000),
        'max_depth': [int(x) for x in np.linspace(10, 110, num=11)] + [None],
        'min_samples_split': randint(2, 50),
        'min_samples_leaf': randint(1, 16),
        'max_features': ['sqrt', 'log2', None, 0.5, 0.8],
        'bootstrap': [True, False],
        'class_weight': [None, 'balanced', 'balanced_subsample'],
        'criterion': ['gini', 'entropy', 'log_loss']
    }

    cv = StratifiedKFold(n_splits=10, shuffle=True, random_state=42)
    model = RandomForestClassifier(random_state=42)  # 初始化 RandomForestClassifier (无需设置 feature_names)
    random_search = RandomizedSearchCV(
        model, param_distributions=param_distributions, n_iter=200, scoring='roc_auc', cv=cv, random_state=42, n_jobs=4
    )
    random_search.fit(X_train_selected, y_train)
    print(f"最佳参数: {random_search.best_params_}")
    print(f"最佳 AUC 得分: {random_search.best_score_}")
    best_model = random_search.best_estimator_

    # 模型评估
    y_pred = best_model.predict(X_val_selected)
    y_pred_proba = best_model.predict_proba(X_val_selected)[:, 1]
    accuracy = accuracy_score(y_val, y_pred)
    print(f"模型准确率: {accuracy}")
    print(classification_report(y_val, y_pred))

    # 绘制并保存 AUC 图 (训练集)
    fpr, tpr, thresholds = roc_curve(y_val, y_pred_proba)
    roc_auc = auc(fpr, tpr)
    plt.figure()
    plt.plot(fpr, tpr, color='darkorange', lw=2, label=f'ROC 曲线 (AUC = {roc_auc:.2f})')
    plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('假阳性率')
    plt.ylabel('真阳性率')
    plt.title('ROC 曲线')
    plt.legend(loc="lower right")
    plt.savefig(f"auc_curve_RF_{method_name}.pdf")  # 保存 AUC 图
    plt.close()

    # 打印验证集的 AUC 值 (也就是绘图展示的 AUC)
    print(f"验证集 AUC 得分: {roc_auc:.2f}")
    # 保存最优模型
    joblib.dump(best_model, f'best_rf_model_{method_name}.pkl')
    print(f"使用 {method_name} 方法的最优模型已保存到 best_rf_model_{method_name}.pkl 文件中。")

    # 预测外部数据
    for i in range(1, 16):
        data_folder = f"./外部数据{i}/批次处理"
        info_file = os.path.join(data_folder, "information.xlsx")
        other_file = os.path.join(data_folder, "other.xlsx")

        # 读取外部数据
        df_info = pd.read_excel(info_file)
        df_other = pd.read_excel(other_file)

        # 预测 other.xlsx 的样本类型
        X_other = df_other.drop(['ID'], axis=1)
        X_other = X_other[train_cols]
        X_other_selected = method.transform(X_other)  # 使用选择的特征进行预测
        df_other['Prediction'] = best_model.predict(X_other_selected)

        # 8. 合并预测结果和生存信息
        df_combined = pd.merge(df_info, df_other[['ID', 'Prediction']], on='ID', how='left')
        df_combined['OS_STATUS'] = df_combined['OS_STATUS'].map({'LIVING': 0, 'DECEASED': 1})

        # 处理缺失值：删除包含缺失值的行
        df_combined.dropna(subset=['OS_MONTHS', 'OS_STATUS', 'Prediction'], inplace=True)

        # 9. 绘制生存曲线并标注 p 值和 HR 值
        # 检查是否有足够的数据绘制生存曲线
        if len(df_combined[df_combined['Prediction'] == 1]) > 0 and len(df_combined[df_combined['Prediction'] == 0]) > 0:
            kmf = KaplanMeierFitter()
            plt.figure()
            groups = ['Non-cluster3', 'Cluster3']
            for k, group in enumerate([0, 1]):
                mask = df_combined['Prediction'] == group
                kmf.fit(df_combined['OS_MONTHS'][mask], df_combined['OS_STATUS'][mask], label=groups[k])
                kmf.plot(ci_show=False)

            results = logrank_test(df_combined[df_combined['Prediction'] == 1]['OS_MONTHS'],
                                   df_combined[df_combined['Prediction'] == 0]['OS_MONTHS'],
                                   event_observed_A=df_combined[df_combined['Prediction'] == 1]['OS_STATUS'],
                                   event_observed_B=df_combined[df_combined['Prediction'] == 0]['OS_STATUS'])
            p_value = results.p_value

            cph = CoxPHFitter(penalizer=0.1)
            cph.fit(df_combined[['OS_MONTHS', 'OS_STATUS', 'Prediction']], duration_col='OS_MONTHS',
                     event_col='OS_STATUS')
            hr = cph.hazard_ratios_['Prediction']

            plt.text(0.1, 0.1, f"HR: {hr:.3f}\np-value: {p_value:.3f}", transform=plt.gca().transAxes)

            plt.title('生存曲线')
            plt.xlabel('时间 (月)')
            plt.ylabel('生存概率')
            plt.savefig(os.path.join(data_folder, f"survival_curve_RF_{method_name}.pdf"))  # 保存生存曲线图
            plt.close()
        else:
            print(f"外部数据{i} 没有足够的样本绘制生存曲线")

        # 10. 保存 df_other 的预测结果
        df_other.to_csv(os.path.join(data_folder, f"other_with_predictions_RF_{method_name}.csv"), index=False)
        print(f"外部数据{i}/other.xlsx 预测结果已保存到 other_with_predictions_RF_{method_name}.csv 文件中。")
    print("-" * 50)  # 分隔不同特征选择方法的结果输出