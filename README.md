# The Data Incubator Semifinal challenge


This repository contains the project proposal for The Data Incubator Fellowship (Semifinal challenge).

## Title: Integration of MRI and genomics data in diagnosing brain aging and Alzheimer’s Disease

### Objectives
1.	To integrate longitudinal brain imaging data and genomics data to identify key features responsible for brain aging and disease
2.	To build a predictive model for diagnosis of Alzheimer’s disease (AD) using machine learning algorithms
3.	To introduce a tool for drug discovery teams to predict response to interventions in cognitive diseases 

### Problem
The 2015 World Alzheimer’s Report suggests that ~50 million people are currently living with dementia, and the rate of growth in dementia cases is exponential. 131.5 million people are projected to suffer from dementia by 2050. Currently, there are 507 ongoing clinical trials to study interventions against dementia and AD, many of which are testing a preventive strategy against these diseases. However, most pharmaceutical companies and clinical investigators lack a robust diagnostic predictive tool, to identify the risk of these diseases and the response to their interventions. Furthermore, individuals who are genetically predisposed to AD or dementia do not have a way to predict their susceptibility to the disease. ![Image](https://github.com/ameya225/TDI/blob/master/Figures/incidence-of-alzheimers.jpg)

### Solution
Fortunately, publicly-accessible large scale datasets with multiple levels of health data are being generated that allow us to utilize sophisticated and novel tools to build predictive models for dementia and AD. Moreover, integrating across different types of data, can not only improve the model accuracy, but also take into account their relationship with each other, since most biological indicators work in unison towards a healthy or a diseased phenotype. Recent advancements in aging research suggest a clear distinction between biological age and chronological age (based on birth year). Hence, regressing variables responsible in age-related diseases like AD/dementia on chronological age is not recommended. A way to address this issue is to use another level of biological information as a surrogate to predict the risk of disease. This surrogate variable can be used to regress the biological variables to build a more accurate predictive model of aging and disease. 

### Data and tools
For my proposed project, I propose to create a surrogate variable dementia using brain imaging data (MRI) obtained from the Open Access Series of Imaging Studies - [OASIS](https://www.oasis-brains.org/). This data is available upon request, but all the images used to build the model can be accessed via my personal GitHub page for this project. To build the predictive model for AD and brain aging, I also plan to use the genetic and immunohistochemistry data from the Aging, Dementia and TBI dataset from the [Allen Brain Atlas](http://aging.brain-map.org/). The normalized gene expression and immunohistochemistry data and metadata is publicly available. Image classification for MRI scans between demented and non-demented donors will be performed using Convolution Neural Networks (CNN) as a part of the Keras-Python library or R Package. The attributes from gene expression and immunohistochemistry data will be selected using the Recursive Feature Elimination tool from the caret-R package.

### Analysis pipeline

#### Image classification steps:
1. The MRI scans from the AD and healthy cohorts are downloaded from OASIS in the NIfTI-1 format. The metadata associated with the participants of the imaging study is downloaded.
2. The NIfTI-1 images are converted into a .jpg format. This is done using [med2image](https://github.com/FNNDSC/med2image) Python3 utility on a high-performance computing cluster (See script). 
3. The .jpg images that 256x256 px will be scaled down to a smaller size (50x50) grayscale images.
4. The image matrix will be then flattened into a vector of length 10000.
5. The dataset will be divided into 75% training and 25% testing partitions.
6. The Keras algorithm will be run- a) Specify architecture; b) Select loss function, optimizer and metrics for goodness-of-fit c) Specify number of epochs, batch size and validation split.
5. The model will be tested, and the efficacy will be computed based on probabilities of accurate predictions.  

##### Healthy brain MRI
![Image](https://github.com/ameya225/TDI/blob/master/OASIS/jpg/OAS30001_ses-d243017.jpg)

##### Dementia brain MRI
![Image](https://github.com/ameya225/TDI/blob/master/OASIS/jpg/OAS30344_ses-d237517.jpg)

#### Genetic and Immunohistochemistry steps 
The entire dataset comprises of 50281 genes from 377 specimens. The 377 specimens are obtained from 107 unique donors. The gene x donor matrix consists of FPKM-Normalized values (fragment-per-kilobase-per-million). FPKM normalization is a way to normalize gene expression data by gene length (in kilobases) and sequencing library size (total number of RNA molecules sequenced).
1. The gene expression data and the immunohistochemistry data will be log normalized and z-transformed (scaled using (x-mean)/sd).
2. Dimensional reduction techniques (PCA/t-SNE) will be used to cluster the entire dataset. The data will be filtered for the top features responsible for the dispersion of the data and classifying the samples into demented and non-demented based on clinical annotation.
3. A generalized linear model will be fit to the gene expression and immunohistochemistry data against the “diagnosis- demented or non-demented” variable with gender, age, race and other clinical categorical variables that are deemed appropriate.
4. A list of top differentially expressed genes between demented and non-demented, based on mean differences will be calculated using a likelihood-ratio test.

#### Integration and Feature selection steps: 
1. The vectors obtained from flattened image data will be merged with the genetic and histological data, by randomly down-sampling the 10000 features.
2. Using R-package- Caret, multiple subset of features will be selected using lasso or elasticnet tools.
3. RMSE values for the model with each subset will be calculated and compared.
4. The model with the lowest RMSE value will be validated on the test data again.

Finally, the entire model will be validated in a different set of MRI and genetic data, before embedding the model onto a user-friendly interface.

### References

A prospective study of cognitive function and onset of dementia in cognitively healthy elders.
Rubin, EH, Storandt, M, Miller, JP, Kinscherf, DA, Grant, EA, Morris, JC, Berg, L, 1998. Arch Neurol. 55, 395-401. PMID: 9520014
