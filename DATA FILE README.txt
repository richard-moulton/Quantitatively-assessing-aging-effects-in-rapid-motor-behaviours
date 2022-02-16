1. Quantitatively-assessing-aging(OH).mat

812 trials are recorded. Each array and struct is synchronized so that all the entries at the ith row of an array (or struct) corresponds to the ith trial.

Attribute Name - Units - Description
====================================
dropsScores - n - The number of distractors successfully "dropped" by a participant.
earlyBinAccuracies - % - The average percentage of targets hit by horizontal bin across all participants' early phases. Phases are determined on a by-participant basis.
earlyBinAvoidances - % - The average percentage of distractors avoided by horizontal bin across all participants' early phases. Phases are determined on a by-participant basis.
hitsScores - n - The number of targets successfully hit by a participant.
maxAvgTargetCreationRate - Hz - The maximum rate at which targets were created during the trial. Calculated as a smoothed average over a 5 s window.
metaDataArray -  - The results of Dexterit-E analysis for the trial. See the Dexterit-E 3.9 User's Guide for more information.
OHflag -  - If 1, then the file contains trials for an Object-Hit task. If 0, then the file contains trials for an Object-Hit-and-Avoid task.
overwhelmedBinAccuracies - % - The average percentage of targets hit by horizontal bin across all participants' overwhelmed phases. Phases are determined on a by-participant basis.
overwhelmedBinAvoidances - % - The average percentage of distractors avoided by horizontal bin across all participants' overwhelmed phases. Phases are determined on a by-participant basis.
steadyStateFrames - # - The frame number at which the participant reached their steady state rate.
steadyState Rates - Hz - The participant's steady state rate for the trial.
participantIDs - # - A participant's unique numerical identifier.
taskScores - 1 - A participant's task score for the trial, as calculated by Dexterit-E. The score itself is a Z-score based off of the results recorded in metaDataArray. See the Dexterit-E 3.9 User's Guide for more information.
trialNames - YYYY-MM-DD_HH-MM-SS - Trials are named by a datetime identifier.
turboFlag = If 1, then the file contains trials for a Turbo task variant. If 0, then the file contains trials for a Classical task variant.

2. Quantitatively-assessing-aging(OH).mat

683 trials are recorded. Each array and struct is synchronized so that all the entries at the ith row of an array (or struct) corresponds to the ith trial.

Attribute Name - Units - Description
====================================
dropsScores - n - The number of distractors successfully "dropped" by a participant.
earlyBinAccuracies - % - The average percentage of targets hit by horizontal bin across all participants' early phases. Phases are determined on a by-participant basis.
earlyBinAvoidances - % - The average percentage of distractors avoided by horizontal bin across all participants' early phases. Phases are determined on a by-participant basis.
hitsScores - n - The number of targets successfully hit by a participant.
maxAvgTargetCreationRate - Hz - The maximum rate at which targets were created during the trial. Calculated as a smoothed average over a 5 s window.
metaDataArray -  - The results of Dexterit-E analysis for the trial. See the Dexterit-E 3.9 User's Guide for more information.
OHflag -  - If 1, then the file contains trials for an Object-Hit task. If 0, then the file contains trials for an Object-Hit-and-Avoid task.
overwhelmedBinAccuracies - % - The average percentage of targets hit by horizontal bin across all participants' overwhelmed phases. Phases are determined on a by-participant basis.
overwhelmedBinAvoidances - % - The average percentage of distractors avoided by horizontal bin across all participants' overwhelmed phases. Phases are determined on a by-participant basis.
steadyStateFrames - # - The frame number at which the participant reached their steady state rate.
steadyState Rates - Hz - The participant's steady state rate for the trial.
participantIDs - # - A participant's unique numerical identifier.
taskScores - 1 - A participant's task score for the trial, as calculated by Dexterit-E. The score itself is a Z-score based off of the results recorded in metaDataArray. See the Dexterit-E 3.9 User's Guide for more information.
trialNames - YYYY-MM-DD_HH-MM-SS - Trials are named by a datetime identifier.
turboFlag = If 1, then the file contains trials for a Turbo task variant. If 0, then the file contains trials for a Classical task variant.

3. Quantitatively-assessing-aging_Participant-demographics.mat

The participants map contains 643 participants mapped to their demographic data using their participant ID as a key. This demographic data consists of 8 fields.
- Sex (coded as M for male or F for female)
- Date of Birth (coded YYYY-MM-DD)
- Handedness (coded as -10 for left-handed, 0 for ambidextrous, and 10 for right handed)
- Age (calculated based on the date of the trial)
- Performed OH (binary flag, 1 indicates performed at least one OH trial)
- Performed OHA (binary flag, 1 indicates performed at least one OHA trial)
- Performed TOH (binary flag, 1 indicates performed at least one TOH trial)
- Performed TOHA (binary flag, 1 indicates performed at least one TOHA trial)

The remaining values are summary statistics for the participants map.

Attribute Name - Units - Description
====================================
maximumAge - yr - The maximum amongst all participant ages.
medianAge - yr - The median of all participant ages.
minimumAge - yr - The minimum amongst all participant ages.
numAmbridextrous - n - The number of participants who are ambidextrous.
numFemale - n - The number of participants who are female.
numLeftHand - n - The number of participants who are left-handed.
numMale - n - The number of participants who are male.
numRightHand - n - The number of participants who are right-handed.
participantsWithErrors - n - A diagnostic counter. If 0, then no errors have been detected in the map.
totalparticipants - n - The total number of participants included in the map.
under16 - n - The total number of participants in the map under the age of 16 years old.
