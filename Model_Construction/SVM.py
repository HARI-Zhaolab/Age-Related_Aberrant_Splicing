import pandas as pd
from sklearn.model_selection import train_test_split, GridSearchCV, StratifiedKFold
from sklearn.svm import SVC
from sklearn.metrics import accuracy_score, classification_report, roc_curve, auc, confusion_matrix, roc_auc_score, f1_score, recall_score, precision_score, matthews_corrcoef
from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import logrank_test
import matplotlib.pyplot as plt
import numpy as np
from sklearn.feature_selection import VarianceThreshold, SelectKBest, mutual_info_classif
from sklearn.svm import LinearSVC
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
import joblib
import os
from sklearn.feature_selection import SelectKBest, f_classif, RFE, SelectFromModel

# 设置中文字体
plt.rcParams['font.sans-serif'] = ['SimHei']
plt.rcParams['axes.unicode_minus'] = False

# 1. 读取训练数据
df_train = pd.read_excel("test.xlsx")

# 2. 数据预处理
def preprocess_data(df):
    # 将 'Cluster3' 映射为 1，其他映射为 0
    df['Group'] = df['Group'].apply(lambda x: 1 if x == 'Cluster3' else 0)
    # 分离特征和目标变量
    X = df.drop(['ID', 'Group'], axis=1)
    y = df['Group']
    return X, y

X, y = preprocess_data(df_train)

# 3. 划分训练集和验证集 (7:3)
X_train, X_val, y_train, y_val = train_test_split(X, y, test_size=0.3, random_state=42)

# 记录训练数据的列名
train_cols = X_train.columns

# 4. 定义特征选择方法，所有方法都选择 10 个特征
feature_selection_methods = {
    #'SelectKBest': SelectKBest(f_classif, k=10),  # 选择基于 F 检验的 10 个最佳特征
    'RFE_LR': RFE(estimator=LogisticRegression(max_iter=100), n_features_to_select=10),  # 使用逻辑回归进行递归特征消除
    'SelectFromModel_RF': SelectFromModel(RandomForestClassifier(n_estimators=100, random_state=42), max_features=10), # 基于随机森林模型的重要性进行特征选择
    'SelectFromModel_LR': SelectFromModel(LogisticRegression(max_iter=100),max_features=10), # 基于逻辑回归模型的 L1 正则化进行特征选择
    'RFE_RF': RFE(estimator=RandomForestClassifier(n_estimators=100, random_state=42), n_features_to_select=10),  # 使用随机森林进行递归特征消除
}

# 5. 训练和评估模型，并预测外部数据
for method_name, method in feature_selection_methods.items():
    print(f"使用特征选择方法: {method_name}")

    # 对训练数据进行特征选择
    X_train_selected = method.fit_transform(X_train, y_train)
    X_val_selected = method.transform(X_val)
    selected_feature_names = X_train.columns[method.get_support()]
    print(f"选择的特征: {selected_feature_names}")

    # SVM 模型参数优化和交叉验证
    param_grid = {
        #'C': np.logspace(-5, 5, 11), # 探索不同的惩罚参数 C
        #'kernel': ['linear', 'poly', 'rbf', 'sigmoid'], # 探索不同的核函数
        #'gamma': ['scale', 'auto'] + list(np.logspace(-7, -1, 7)),# 探索不同的 gamma 参数
        'class_weight': [None, 'balanced']# 探索不同的类别权重
    }

    cv = StratifiedKFold(n_splits=10, shuffle=True, random_state=42) # 10 折分层交叉验证
    model = SVC(probability=True, random_state=42) # 初始化 SVM 模型
    grid_search = GridSearchCV(model, param_grid, scoring='roc_auc', cv=cv, n_jobs= 40) # 使用 GridSearchCV 进行参数优化
    grid_search.fit(X_train_selected, y_train) # 拟合模型

    print(f"最佳参数: {grid_search.best_params_}")
    print(f"最佳 AUC 得分: {grid_search.best_score_}")
    best_model = grid_search.best_estimator_ # 获取最佳模型

    # 模型评估
    y_pred = best_model.predict(X_val_selected) # 预测验证集标签
    y_pred_proba = best_model.predict_proba(X_val_selected)[:, 1] # 预测验证集概率
    accuracy = accuracy_score(y_val, y_pred) # 计算准确率

    # 计算混淆矩阵
    cm = confusion_matrix(y_val, y_pred)

    # 计算其他指标
    auc_score = roc_auc_score(y_val, y_pred_proba)
    f1 = f1_score(y_val, y_pred)
    recall = recall_score(y_val, y_pred)
    precision = precision_score(y_val, y_pred)
    mcc = matthews_corrcoef(y_val, y_pred)  # 计算 MCC

    # 打印验证集的指标
    print(f"验证集指标：")
    print(f"模型准确率: {accuracy}")
    print(f"混淆矩阵:\n{cm}")
    print(f"AUC 得分: {auc_score:.2f}")
    print(f"F1 得分: {f1:.2f}")
    print(f"召回率: {recall:.2f}")
    print(f"精确率: {precision:.2f}")
    print(f"马修斯相关系数 (MCC): {mcc:.2f}") # 打印 MCC
    print(classification_report(y_val, y_pred))

    # 绘制 AUC 图 (训练集)
    fpr, tpr, thresholds = roc_curve(y_val, y_pred_proba) # 计算 ROC 曲线
    roc_auc = auc(fpr, tpr) # 计算 AUC
    plt.figure()
    plt.plot(fpr, tpr, color='darkorange', lw=2, label=f'ROC 曲线 (AUC = {roc_auc:.2f})')
    plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('假阳性率')
    plt.ylabel('真阳性率')
    plt.title('ROC 曲线')
    plt.legend(loc="lower right")
    plt.savefig(f"auc_curve_SVM_{method_name}.pdf") # 保存 AUC 图
    plt.close()

    # 打印验证集的 AUC 值 (也就是绘图展示的 AUC)
    print(f"验证集 AUC 得分: {roc_auc:.2f}")
    # 保存最优模型
    joblib.dump(best_model, f'best_svm_model_{method_name}.pkl')
    print(f"使用 {method_name} 方法的最优模型已保存到 best_svm_model_{method_name}.pkl 文件中。")

    # 预测外部数据
    for i in range(1, 16):
        data_folder = f"./外部数据{i}/批次处理"
        info_file = os.path.join(data_folder, "information.xlsx")
        other_file = os.path.join(data_folder, "other.xlsx")

        # 读取外部数据
        df_info = pd.read_excel(info_file)
        df_other = pd.read_excel(other_file)

        # 预测 other.xlsx 的样本类型
        X_other = df_other.drop(['ID'], axis=1) # 删除 'ID' 列

        # 使用训练数据列名确保特征顺序一致
        X_other = X_other[train_cols]

        X_other_selected = method.transform(X_other)  # 使用选择的特征进行预测
        df_other['Prediction'] = best_model.predict(X_other_selected) # 预测外部数据标签

        # 9. 合并预测结果和生存信息
        df_combined = pd.merge(df_info, df_other[['ID', 'Prediction']], on='ID', how='left')
        df_combined['OS_STATUS'] = df_combined['OS_STATUS'].map({'LIVING': 0, 'DECEASED': 1}) # 将生存状态映射为 0 和 1

        # 处理缺失值：删除包含缺失值的行
        df_combined.dropna(subset=['OS_MONTHS', 'OS_STATUS', 'Prediction'], inplace=True)

        # 10. 绘制生存曲线并标注 p 值和 HR 值
        # 检查是否有足够的数据绘制生存曲线
        group_0_count = df_combined[df_combined['Prediction'] == 0]['OS_MONTHS'].count() # 计算预测为 0 的样本数量
        group_1_count = df_combined[df_combined['Prediction'] == 1]['OS_MONTHS'].count() # 计算预测为 1 的样本数量

        if group_0_count >= 5 and group_1_count >= 5: # 确保每组至少有 5 个样本
            kmf = KaplanMeierFitter()
            plt.figure()
            groups = ['Non-cluster3', 'Cluster3']
            for k, group in enumerate([0, 1]):
                mask = df_combined['Prediction'] == group
                kmf.fit(df_combined['OS_MONTHS'][mask], df_combined['OS_STATUS'][mask], label=groups[k]) # 拟合 Kaplan-Meier 模型
                kmf.plot(ci_show=False) # 绘制生存曲线

            results = logrank_test(df_combined[df_combined['Prediction'] == 1]['OS_MONTHS'],
                                   df_combined[df_combined['Prediction'] == 0]['OS_MONTHS'],
                                   event_observed_A=df_combined[df_combined['Prediction'] == 1]['OS_STATUS'],
                                   event_observed_B=df_combined[df_combined['Prediction'] == 0]['OS_STATUS']) # 进行对数秩检验
            p_value = results.p_value # 获取 p 值

            cph = CoxPHFitter(penalizer=0.1) # 初始化 CoxPH 模型
            cph.fit(df_combined[['OS_MONTHS', 'OS_STATUS', 'Prediction']], duration_col='OS_MONTHS',
                    event_col='OS_STATUS') # 拟合 COX 模型
            hr = cph.hazard_ratios_['Prediction'] # 获取风险比

            plt.text(0.1, 0.1, f"HR: {hr:.3f}\np-value: {p_value:.3f}", transform=plt.gca().transAxes) # 在图上添加风险比和 p 值

            plt.title('生存曲线')
            plt.xlabel('时间 (月)')
            plt.ylabel('生存概率')
            plt.savefig(f"./外部数据{i}/批次处理/survival_curve_SVM_{method_name}.pdf") # 保存生存曲线图
            plt.close()
        else:
            print(
                f"外部数据 {i} 的分组样本数不足 (Non-cluster3: {group_0_count}, Cluster3: {group_1_count})，无法绘制生存曲线。")

        # 11. 保存 df_other 的预测结果
        df_other.to_excel(f"./外部数据{i}/批次处理/other_with_predictions_SVM_{method_name}.xlsx", index=False)
        print(
            f"df_other 的预测结果已保存到 ./外部数据{i}/批次处理/other_with_predictions_SVM_{method_name}.xlsx 文件中。")
    print("-" * 50)