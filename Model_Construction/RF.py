import pandas as pd
from sklearn.model_selection import train_test_split, RandomizedSearchCV, StratifiedKFold
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score, classification_report, roc_curve, auc
import matplotlib.pyplot as plt
import numpy as np

# 1. 读取数据
df = pd.read_csv("text.csv")

df_info = pd.read_excel("information.xlsx")
df_other = pd.read_excel("other.xlsx")

# 2. 数据预处理
df['Group'] = df['Group'].apply(lambda x: 1 if x == 'Cluster3' else 0)
X = df.drop(['ID', 'Group'], axis=1)
y = df['Group']

# 3. 划分训练集和测试集 (7:3)
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=42)

# 4. 随机森林模型参数优化和交叉验证
# 定义更详细的参数空间，使用随机搜索
param_distributions = {
    'n_estimators': [int(x) for x in np.linspace(start=100, stop=1000, num=10)],
    'max_depth': [int(x) for x in np.linspace(10, 110, num=11)] + [None],
    'min_samples_split': [2, 5, 10, 20, 50],
    'min_samples_leaf': [1, 2, 4, 8, 16],
    'max_features': ['sqrt', 'log2', None, 0.5, 0.8],  # 添加更多max_features选项
    'bootstrap': [True, False],
    'class_weight': [None, 'balanced', 'balanced_subsample'],
    'criterion': ['gini', 'entropy', 'log_loss']  # 添加criterion参数
}

# 使用 StratifiedKFold 进行分层交叉验证，增加折数
cv = StratifiedKFold(n_splits=10, shuffle=True, random_state=42)

# 使用 RandomizedSearchCV 进行随机搜索，增加迭代次数
model = RandomForestClassifier(random_state=42)
random_search = RandomizedSearchCV(
    model, param_distributions=param_distributions, n_iter=200, scoring='roc_auc', cv=cv, random_state=42, n_jobs=-1
)
random_search.fit(X_train, y_train)

# 输出最佳参数和最佳得分
print(f"最佳参数: {random_search.best_params_}")
print(f"最佳 AUC 得分: {random_search.best_score_}")

# 使用最佳模型进行预测
best_model = random_search.best_estimator_
y_pred = best_model.predict(X_test)
y_pred_proba = best_model.predict_proba(X_test)[:, 1]

# 5. 模型评估
accuracy = accuracy_score(y_test, y_pred)
print(f"模型准确率: {accuracy}")
print(classification_report(y_test, y_pred))

nan_predictions = df_other[df_other['Prediction'].isnull()]
print(nan_predictions)
# 6. 绘制 AUC 图
fpr, tpr, thresholds = roc_curve(y_test, y_pred_proba)
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
plt.show()

# 7. 预测所有样本的分组
all_predictions = best_model.predict(X)

# 8. 将预测结果添加到原数据中
df['Prediction'] = all_predictions

# 9. 保存预测结果
df.to_csv("text_with_predictions_optimized.csv", index=False)

print("预测结果已保存到 text_with_predictions_optimized.csv 文件中。")