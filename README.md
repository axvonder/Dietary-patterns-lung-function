Hi. Thanks for accessing the files of our project.

Project title: Is dietary pattern associated with rate of decline in lung function?

Project Description

This project investigates the association between dietary patterns and lung function decline, focusing on the Healthy Eating Index-2020 (HEI-2020) and the Dietary Approaches to Stop Hypertension (DASH) index. We analyzed longitudinal data from two large cohort studies to understand how diet impacts lung function parameters, such as forced expiratory volume (FEV1) and forced vital capacity (FVC), and their decline over time. The study also explores the effects of smoking status and demographic factors on these associations. Our findings contribute to knowledge about dietary interventions for maintaining lung health and preventing chronic conditions like COPD.

Table of Contents

1. Introduction
2. Data Source
3. Repository Structure
4. Installation and Setup
5. Usage
6. Contributing
7. License
8. Contact
9. Acknowledgments

## 1. Introduction

Authors: Alexander Vonderschmidt, Bonnie K Patchen, Kathryn Bass, Kathryn B Arnold, Eric J Shiroma, Eleanor Simonsick, Dana B Hancock, and Patricia A Cassano. The publication can be found [HERE-CITE PUBLICATION WHEN PUBLIC].
The main objectives of the analysis are:

- Investigate associations between dietary patterns (HEI-2020 and DASH) and longitudinal change in lung function, as measured by: FEV1, FVC, and FEV1/FVC ratio.
- Explore how these changes differ by smoking status and key demographic characteristics such as sex and race.

## 2. Data Source

The analysis uses data from the Health, Aging, and Body Composition Study (Health ABC) and the Respiratory Ancillary Study (RAS), a sub-study of the Selenium and Vitamin E Cancer Prevention Trial (SELECT). The Health ABC data is available upon request through the NIH's NIA. To inquire about data sharing agreements with RAS, please contact the corresponding author, PAC.

## 3. Repository Structure

README.md: This file, providing an overview of the project as well as instructions on installation and usage.
FFQ coding.xlsx: A multi-sheet excel workbook containing the breakdown of the food frequency questionnaire (FFQ) items into the HEI-2020 and DASH categories for each of the two cohorts. For each FFQ item, if the serving/portion size was not directly listed within the FFQ, we used USDA databases to source an estimated portion size. Each estimation has relevant links directly in-line for each FFQ in which this applies. Additionally, some composite foods (such as “stew, pot pie, curries and casseroles with meat or chicken” in the RAS FFQ) were labeled internally as 'mixed dishes' and a percent distribution to many food categories was applies. All mixed dish estimations are also listed in-line for each FFQ item in which this applies.
Dietary pattern indices code (HEI2020, DASH).R: An R script that takes in the raw data from Health ABC and RAS, cleans and processes the data to create the variables needed for analysis, and applies the aforementioned FFQ coding into the relevant dietary pattern scores.
Dietary patterns analyses UPDATE July 2024.R: This file handles all of the analyses presented in the manuscript. It will take the data files—in .csv format—output from "Dietary pattern indices code (HEI2020, DASH).R" and output 1 excel workbook containing all tables in the manuscript in individual sheets within the workbook. It will also output the figures presented in the manuscript as .png files.

## 4. Installation and Setup

1. Clone this repository to your local machine:
git clone https://github.com/axvonder/Dietary-patterns-lung-function.git
Install the required R packages:
install.packages(c("lme4", "metafor", "survey", "tidyverse", "ggplot2", "srvyr"))
Set your working directory to the cloned repository folder in your R environment.

## 5. Usage

Run the R scripts in the following order:

1. Dietary pattern indices code (HEI2020, DASH).R: cleans, processes, transforms data for analysis, outputs .csv files to use for analysis.
2. Dietary patterns analyses cornell.R: runs all analyses, produces 1 excel workbook with each table as a separate sheet. Produces the figures used in the manuscript as .png files.

## 6. Contributing

Contributions to this project are welcome. If you would like to contribute, please follow these steps:
1. Fork the repository to your GitHub account.
2. Clone your fork to your local machine: `https://github.com/axvonder/NDNSMeatTrends.git`
3. Create a new branch for your changes: `git checkout -b feature/my-new-feature`
4. Make your changes and commit them to your branch: `git add .` and `git commit -m "Add my new feature"`
5. Push your changes to your fork: `git push origin feature/my-new-feature`
6. Open a pull request on the original repository with a clear and concise description of your changes.

Please ensure that your changes are consistent with the project's style and that you have tested your code before submitting a pull request. Also, include any relevant documentation or comments in your code.

## 8. Contact

For any questions or concerns, please contact the senior and corresponding author, PAC. Alternatively, you can contact the project maintainer:

- Alexander Vonderschmidt: a.vonderschmidt@sms.ed.ac.uk

## 9. Acknowledgments

We thank all co-authors and collaborators, including the participants and investigators of the Health ABC and RAS cohorts.
