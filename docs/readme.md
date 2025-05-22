# panc8_scina Project Documentation

## Project Structure

<img width="219" alt="image" src="https://github.com/user-attachments/assets/a9f418ad-9bf5-4780-bd29-b70d4fe40884" />



## Project Directory Explanation

### docs
This directory is used to store project documentation, such as project description documents and user manuals.

### input
This directory stores the input data for the project, currently containing the `panc8_scina.h5ad` file, which is a single-cell RNA sequencing data file and the main input data for the project analysis.

### models
This directory is used to store model files used in the project, such as `scimilarity` model files, which are used for cell annotation tasks in the project.

### notebooks
This directory stores Jupyter Notebook files, currently containing `panc8_scina.ipynb`, which is the main analysis code file of the project. It records a series of analysis processes from data loading, preprocessing, cell annotation to result visualization.

### output
This directory stores the output results of the project, divided into two subdirectories:

- `output_modified_h5ad`: Stores modified h5ad files, such as `panc8_scina_mod.h5ad`, which are generated after processing and annotating the original input data. These files contain new annotation information.
- `output_report`: Stores analysis result report files, such as `result_panc8_scina.csv`, which records evaluation metrics (such as precision, recall, F1 score) of cell annotation results.

### scripts
This directory stores project-related script files, currently containing `download_model.sh`, which may be used to download model files and other resources required by the project.


## panc8_scina.ipynb File Description

### Data Loading and Preprocessing
- Use the `scanpy` library to load the single-cell RNA sequencing data file `panc8_scina.h5ad` and set relevant plotting parameters and ignore warning messages.
- Import the `CellAnnotation` class and related utility functions from the `scimilarity` library for subsequent cell annotation tasks.
- Align the dataset with the feature space of the `scimilarity` model to ensure data consistency.
- Normalize the data (the call to the `lognorm_counts` function is commented out and can be enabled as needed).

### Cell Annotation
- Use the `scimilarity` model to calculate the embedding vectors of the data and visualize them using the UMAP algorithm to show the distribution of cells.
- Perform unconstrained annotation, classifying cells into any type present in the `scimilarity` reference, and visualize the annotation results. Filter out well-represented cell types in the data for display.
- Conduct constrained classification by specifying a list of target cell types, use the `scimilarity` model to annotate the dataset, and visualize the annotation results.

### Annotation Quality Assessment
- Visualize the quality of the annotation results using UMAP by displaying the distance of each cell to the nearest reference cell in the `min_dist` column to evaluate annotation accuracy.
- Extract the actual cell types and predicted cell types for comparison and analysis. Count true positives (TP), false positives (FP), and false negatives (FN), and calculate evaluation metrics such as precision, recall, and F1 score. Save the results to the `result_panc8_scina.csv` file.

### Output Results
- Save the processed and annotated data as a new h5ad file `panc8_scina_mod.h5ad` for further analysis and use.

## Usage Instructions

1. **Ensure Project Structure**: Make sure the project directory structure is complete and all related files and data are prepared.
2. **Install Dependencies**: Install the required Python dependencies for the project, such as `scanpy` and `scimilarity`.
3. **Run the Notebook**:
   - Navigate to the `notebooks` directory.
   - Open and run the `panc8_scina.ipynb` file. Execute the cells in order to complete data loading, preprocessing, cell annotation, quality assessment, and result output.
4. **Check Outputs**:
   - Navigate to the `output` directory.
   - Check the output files, such as `panc8_scina_mod.h5ad` and `result_panc8_scina.csv`, to obtain the analysis results and evaluation metrics.
  

## Notes

- **Path Configuration**: Before running the code, ensure that the paths to the data file `panc8_scina.h5ad` and model files are correct. Modify the path configurations in the code if necessary.
- **Extending Functionality**: If further analysis of the data or extension of functionality is required, additional code cells can be added to the `panc8_scina.ipynb` file for implementation.


## Machine Requirements
1. **RAM**: For Windows and Linux users, 32GB and more are required. For Mac users, 16GB is enough.
2. **ROM**: At least 100GB of free space required.
3. **Operating System**: Unix-like systems are preferred, since some features are hard encoded.


