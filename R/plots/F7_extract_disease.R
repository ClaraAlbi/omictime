

library(tibble)

disease_icd10 <- tribble(
  ~condition, ~vocab_id, ~code, ~code_desc,

  # --- Neuro-psychiatric ---
  "rec_MDD_mild", "ICD10", "F33.0", "Recurrent depressive disorder, current episode mild",
  "rec_MDD_moderate", "ICD10", "F33.1", "Recurrent depressive disorder, current episode moderate",
  "rec_MDD_severe", "ICD10", "F33.2", "Recurrent depressive disorder, current episode severe without psychotic symptoms",
  "rec_MDD_severe_with_psychosis", "ICD10", "F33.3", "Recurrent depressive disorder, current episode severe with psychotic symptoms",
  "rec_MDD_remission", "ICD10", "F33.4", "Recurrent depressive disorder, currently in remission",
  "rec_MDD_other", "ICD10", "F33.8", "Other recurrent depressive disorders",
  "rec_MDD_unspecified", "ICD10", "F33.9", "Recurrent depressive disorder, unspecified",

  "episode_MDD_mild", "ICD10", "F32.0", "Mild depressive episode",
  "episode_MDD_moderate", "ICD10", "F32.1", "Moderate depressive episode",
  "episode_MDD_severe", "ICD10", "F32.2", "Severe depressive episode without psychotic symptoms",
  "episode_MDD_severe_with_psychosis", "ICD10", "F32.3", "Severe depressive episode with psychotic symptoms",
  "episode_MDD_other", "ICD10", "F32.8", "Other depressive episodes",
  "episode_MDD_unspecified", "ICD10", "F32.9", "Depressive episode, unspecified")

"bipolar_disorder", "ICD10", "F31.0", "Bipolar affective disorder, current episode hypomanic",
"bipolar_disorder", "ICD10", "F31.1", "Bipolar affective disorder, current episode manic without psychotic symptoms",
"bipolar_disorder", "ICD10", "F31.2", "Bipolar affective disorder, current episode manic with psychotic symptoms",
"bipolar_disorder", "ICD10", "F31.3", "Bipolar affective disorder, current episode mild or moderate depression",
"bipolar_disorder", "ICD10", "F31.4", "Bipolar affective disorder, current episode severe depression without psychotic symptoms",
"bipolar_disorder", "ICD10", "F31.5", "Bipolar affective disorder, current episode severe depression with psychotic symptoms",
"bipolar_disorder", "ICD10", "F31.6", "Bipolar affective disorder, currently in remission",
"bipolar_disorder", "ICD10", "F31.7", "Bipolar affective disorder, currently unspecified",
"bipolar_disorder", "ICD10", "F31.8", "Other bipolar affective disorders",
"bipolar_disorder", "ICD10", "F31.9", "Bipolar affective disorder, unspecified",

"anxiety", "ICD10", "F41.0", "Panic disorder [episodic paroxysmal anxiety]",
"anxiety", "ICD10", "F41.1", "Generalized anxiety disorder",
"anxiety", "ICD10", "F41.2", "Mixed anxiety and depressive disorder",
"anxiety", "ICD10", "F41.3", "Other mixed anxiety disorders",
"anxiety", "ICD10", "F41.8", "Other specified anxiety disorders",
"anxiety", "ICD10", "F41.9", "Anxiety disorder, unspecified",

# --- Cardiometabolic ---
"type2_diabetes", "ICD10", "E11.0", "Type 2 diabetes mellitus with coma",
"type2_diabetes", "ICD10", "E11.1", "Type 2 diabetes mellitus with ketoacidosis",
"type2_diabetes", "ICD10", "E11.2", "Type 2 diabetes mellitus with kidney complications",
"type2_diabetes", "ICD10", "E11.3", "Type 2 diabetes mellitus with ophthalmic complications",
"type2_diabetes", "ICD10", "E11.4", "Type 2 diabetes mellitus with neurological complications",
"type2_diabetes", "ICD10", "E11.5", "Type 2 diabetes mellitus with circulatory complications",
"type2_diabetes", "ICD10", "E11.6", "Type 2 diabetes mellitus with other specified complications",
"type2_diabetes", "ICD10", "E11.7", "Type 2 diabetes mellitus with multiple complications",
"type2_diabetes", "ICD10", "E11.8", "Type 2 diabetes mellitus with unspecified complications",
"type2_diabetes", "ICD10", "E11.9", "Type 2 diabetes mellitus without complications",

"hypertension", "ICD10", "I10", "Essential (primary) hypertension",
"hypertension", "ICD10", "I11", "Hypertensive heart disease",
"hypertension", "ICD10", "I12", "Hypertensive renal disease",
"hypertension", "ICD10", "I13", "Hypertensive heart and renal disease",
"hypertension", "ICD10", "I15", "Secondary hypertension",

# --- Immune ---
"rheumatoid_arthritis", "ICD10", "M05.0", "Felty's syndrome",
"rheumatoid_arthritis", "ICD10", "M05.1", "Rheumatoid lung disease",
"rheumatoid_arthritis", "ICD10", "M05.2", "Rheumatoid vasculitis",
"rheumatoid_arthritis", "ICD10", "M05.3", "Rheumatoid arthritis with involvement of other organs or systems",
"rheumatoid_arthritis", "ICD10", "M05.8", "Other seropositive rheumatoid arthritis",
"rheumatoid_arthritis", "ICD10", "M05.9", "Seropositive rheumatoid arthritis, unspecified",
"rheumatoid_arthritis", "ICD10", "M06.0", "Seronegative rheumatoid arthritis",
"rheumatoid_arthritis", "ICD10", "M06.9", "Rheumatoid arthritis, unspecified",

"asthma", "ICD10", "J45.0", "Predominantly allergic asthma",
"asthma", "ICD10", "J45.1", "Non-allergic asthma",
"asthma", "ICD10", "J45.8", "Mixed asthma",
"asthma", "ICD10", "J45.9", "Asthma, unspecified",
"asthma", "ICD10", "J46", "Status asthmaticus",

# --- Cancer ---
"breast_cancer", "ICD10", "C50.0", "Malignant neoplasm of nipple and areola",
"breast_cancer", "ICD10", "C50.1", "Malignant neoplasm of central portion of breast",
"breast_cancer", "ICD10", "C50.2", "Malignant neoplasm of upper-inner quadrant of breast",
"breast_cancer", "ICD10", "C50.3", "Malignant neoplasm of lower-inner quadrant of breast",
"breast_cancer", "ICD10", "C50.4", "Malignant neoplasm of upper-outer quadrant of breast",
"breast_cancer", "ICD10", "C50.5", "Malignant neoplasm of lower-outer quadrant of breast",
"breast_cancer", "ICD10", "C50.8", "Malignant neoplasm of overlapping sites of breast",
"breast_cancer", "ICD10", "C50.9", "Malignant neoplasm of breast, unspecified"
)
