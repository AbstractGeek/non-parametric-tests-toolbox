# Nonparametric tests (MATLAB functions)
A toolbox of standard non-parametric tests (MATLAB).

## Mann Whitney U test
Performs the mann Whitney U test on the two groups of data. This is the nonparametric equivalent of the student t-test. Refer to https://en.wikipedia.org/wiki/Mann%E2%80%93Whitney_U_test for more general details.

The mathematical formulation of the test for this function was based on "Biostatistical Analysis" by Jerrold H Zar.

## Kruskal Wallis test
Performs the Kruskal Wallis one-way analysis of variance on all the groups (samples). This is the nonparametric version of one-way anova. Refer to https://en.wikipedia.org/wiki/Kruskal%E2%80%93Wallis_one-way_analysis_of_variance for more general details.

The mathematical formulation of the test for this function was based on Section 11.6 (Nonparametric multiple comparisons) in "Biostatistical Analysis" by Jerrold H Zar.

## Nemenyi test
This is the post-hoc test used to find groups of data that defer from each other. This test is used after Kruskal-Wallis test rejects the null hypothesis. Nemenyi test is a nonparametric version of the Tukey test. Refer to "Biostatistical Analysis" by Jerrold H Zar for more details.

The mathematical formulation of the test for this function was based on Section 11.6 (Nonparametric multiple comparisons) in "Biostatistical Analysis" by Jerrold H Zar.

## PerformStats
PerformStats function performs both Kruskal-Wallis and Nemenyi test on groups of data. It is a convenient wrapper function to quickly perform all pertaining nonparametric tests on groups of data.
