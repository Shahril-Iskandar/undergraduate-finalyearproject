# undergraduate-finalyearproject
This repository contains the scripts used during my undergraduate final year honors thesis project.

Project description:
- Aim: Assess the accuracy of FLEX, a velocity-based training device, in measuring mean and peak concentric velocity in 2 weightlifting movements (power clean and snatch) during the concentric phase.
- Criterion device used: VICON 3D motion capture (mocap)
- Methods:
- Participants performed varying repetitions across the load intensity (20-90% of 1RM) for each exercise.
- Data was captured simultaneously from FLEX and 3D mocap
- Bland-Altman was used to calculate the mean differences and limits of agreement.

ForVICONv3.m
- Analyzed each repetition time-series data from the reflective markers placed on both ends of the barbell to extract the mean and peak concentric velocity during the concentric phase.

ValidityTesting.m
- Perform the Bland-Altman analysis
