#' Covert SIA activities into 4 groupings and exclude (or not) based on logicals
#'
#' @param df_sia - dataframe
#' @param PIRI - logical - to include SIAs Activities focused on RI or not
#' @param OR - logical - to include outbreak response SIAs Activities or not
#' @param OTHER - logical - to include any other SIAs Activities or not
#'
#' @import tidyverse
#'
#' @returns list of new dataframe AND vector of excluded SIAs
ConvertActivitySIA <- function(df_sia, PIRI, OR, OTHER){

  # Broke SIAs into four types (PIRI, OR, OTHER, PREVENTATIVE)
  #PIRI = "SNID", "NID", "NIDs","VACWEEK", "PIRI", "Vaccination Week", "VaccinationWeek", "CHD", "Routine", "Vaccination week"
  #OR = "Emergency campaign", "FollowUp/OR", "OR", "CaseResponse", "Case response", "Outbreak Response", "Case Response"
  #OTHER = "HighRiskAreas", "Defaulter Tracing", "Ring vaccination","UNDEFINED", NA, "HRA", "Other", "MR campaign among adult risk group age 20-40 Y", "High risk Areas", "HRAs", "High Risk Areas"
  #PREVENTATIVE = "Campaign","CatchUp","FollowUp","MopUp","SpeedUp","Followup","FollowUp - Phased","CatchUp-SIA","Preventive","MassPreventive","Mass Preventive","SIA", "Speedup", "FILLIN", "preventive", "Mass campaign", "Catch-up SIA", "FollowUp or CatchUp as per JRF?"

  SIA_include_opts <- c("PREVENTATIVE","PIRI","OR","OTHER")
  SIA_include <- SIA_include_opts[c(TRUE,PIRI,OR,OTHER)]

  # This code makes sure the above are captured correctly, plus hopefully any new version
  df_sia <- df_sia |>
    dplyr::mutate(
      Activity_clean = stringr::str_to_lower(ACTIVITY_TYPE),
      Activity_clean = stringr::str_trim(Activity_clean),
      Activity_Cat = dplyr::case_when(
        is.na(Activity_clean) ~ "OTHER",
        stringr::str_detect(Activity_clean,"nid|snid|week|piri|chd|routine") ~ "PIRI",
        stringr::str_detect(Activity_clean,"hra|defaulter|ring|undefined|risk") ~ "OTHER",
        stringr::str_detect(Activity_clean,"outbreak|case.?response|^or$|emergency|follow.?up/or") ~ "OR",
        stringr::str_detect(Activity_clean,"campaign|catch.?up|follow.?up|mop.?up|speed.?up|prevent|mass|sia|fill") ~ "PREVENTATIVE",
        TRUE ~ "PREVENTATIVE"
      )
    )

  df_sia1 <- df_sia |>
    dplyr::filter((Activity_Cat %in% SIA_include)) |>
    dplyr::select(-Activity_clean, -Activity_Cat)

  sia_pulled <- df_sia |>
    dplyr::filter(!(Activity_Cat %in% SIA_include)) |>
    dplyr::pull(Activity_clean) |> unique()

  return(df_sia1)

}
