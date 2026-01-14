import pandas as pd
import numpy as np
from sklearn.linear_model import LinearRegression, LogisticRegression
from doubleml.data import DoubleMLPanelData, DoubleMLDIDData, DoubleMLData
from doubleml.did import DoubleMLDIDMulti
import os

os.chdir(
    "/Users/laurensmauz/Documents/causal_machine_learning/Group_assignment/data_code_cash_bail_study"
)

dta_study = pd.read_stata("CashBail.dta")
dta_study["EligibleOffense"] = dta_study["EligibleOffense"].astype("float")
dta_study.loc[dta_study["EligibleOffense"] == 0, "EligibleOffense"] = (
    np.inf
)  # this is done so that doubleML package knows the never-treated obs.
dta_study_clean = dta_study[
    [
        "docketnumber_anon",
        "EligibleOffense",
        "Post",
        "hasfta",
        "defendantageatarrest",
        "male",
        "defendantisblack",
        "Hisp",
        "prior",
        "prior_FTA",
        "felony",
        "nb_charges",
    ]
]


dml_data = DoubleMLPanelData(
    data=dta_study_clean,
    y_col="hasfta",
    d_cols="EligibleOffense",
    t_col="Post",
    x_cols=[
        "defendantageatarrest",
        "male",
        "defendantisblack",
        "Hisp",
        "prior",
        "prior_FTA",
        "felony",
        "nb_charges",
    ],
    id_col="docketnumber_anon",
)
print(dml_data)

dml_obj = DoubleMLDIDMulti(
    obj_dml_data=dml_data,
    ml_g=LinearRegression(),
    ml_m=LogisticRegression(max_iter=1000),
    control_group="never_treated",
    panel=False,
)

dml_obj.fit()
print(dml_obj)
