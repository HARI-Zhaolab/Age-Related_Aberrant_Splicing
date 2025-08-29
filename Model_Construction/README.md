# Model_Construction

This project contains scripts for building machine learning classification models and external validation datasets from multiple breast cancer cohorts. The main goal is to demonstrate the prognostic prediction capability of the models.

------

## Project Structure

### Datasets

- **GSE20685/**
- **GSE25066/**
- **GSE42568/**
- **Metabric/**

These folders contain publicly available breast cancer datasets (gene expression / clinical information / prediction results from different models).

------

### Model & Algorithm Scripts

- **DT.py**
   Classification using Decision Tree.
- **DA.py**
   Discriminant analysis methods (e.g., LDA/QDA).
- **GP.py**
   Gaussian Process for classification/regression.
- **XGBoost.py**
   Gradient Boosting Trees model.
- **KNN.py**
   K-Nearest Neighbors classifier.
- **LD.py**
   Linear Discriminant model (LDA).
- **RF.py**
   Random Forest training and evaluation.
- **NB.py**
   Naïve Bayes classifier.
- **LR.py**
   Logistic Regression model.
- **SVM.py**
   Support Vector Machine model.

------

### Saved Models

- **best_svm_model_SelectKBest.pkl**
   Pre-trained **SVM + SelectKBest feature selection** model saved as a pickle file, ready to load and use.

------

### Documentation

- **README.md**
   Project description file.
- **requirements.txt**
   Environment and dependency information for reproducibility.

------

## Usage

1. Prepare the required raw data (place it in the corresponding GSE dataset folder following the provided format).
2. Run the desired script, for example:

```python
bash
python SVM.py
```