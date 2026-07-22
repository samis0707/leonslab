"""Volume calculations for the DQ Elastin elastase assay."""
from __future__ import annotations

# --- Konstanten ---
SAMPLE_VOL_PER = 10
SAMPLE_RB_PER = 590
NEG_CTRL_LB_VOL = 10
NEG_CTRL_RB = 590

ELASTASE_STOCK_VOL = 5
STD_TUBE1_RB = 745
DIL_TRANSFER_VOL = 375
DIL_RB_PER_STEP = 375
NUM_DILUTIONS = 6
STD_ZERO_RB = 375
NUM_STANDARDS = 8

VOL_PER_WELL = 150
DUPLICATES = 2
NUM_NEG_CTRL = 1

DQ_WORKING_PER_WELL = 50
DQ_STOCK_PER_WORKING = 5
DQ_RB_PER_WORKING = 45
DQ_OVERSHOOT = 0.20

RB_RESERVE = 0.10


def calculate_elastase_assay(num_samples: int) -> dict:
    """Berechnet alle benötigten Volumina für einen Elastase-Assay."""
    wells_std = NUM_STANDARDS * DUPLICATES
    wells_samp = num_samples * DUPLICATES
    wells_neg = NUM_NEG_CTRL * DUPLICATES
    wells_total = wells_std + wells_samp + wells_neg

    rb_samples = SAMPLE_RB_PER * num_samples
    rb_neg_ctrl = NEG_CTRL_RB

    rb_std_tube1 = STD_TUBE1_RB
    rb_std_dilutions = DIL_RB_PER_STEP * NUM_DILUTIONS
    rb_std_zero = STD_ZERO_RB
    rb_standard = rb_std_tube1 + rb_std_dilutions + rb_std_zero

    dq_working_needed = wells_total * DQ_WORKING_PER_WELL
    dq_working_total = dq_working_needed * (1 + DQ_OVERSHOOT)
    dq_stock_needed = dq_working_total * (DQ_STOCK_PER_WORKING / DQ_WORKING_PER_WELL)
    rb_dq_working = dq_working_total * (DQ_RB_PER_WORKING / DQ_WORKING_PER_WELL)

    rb_1x_subtotal = rb_samples + rb_neg_ctrl + rb_standard + rb_dq_working
    rb_1x_total = rb_1x_subtotal * (1 + RB_RESERVE)

    rb_10x_needed = rb_1x_total / 10
    h2o_needed = rb_1x_total * 9 / 10

    concentrations = []
    base_conc = 0.5
    for i in range(7):
        concentrations.append(base_conc / (2**i))
    concentrations.append(0.0)

    return {
        "num_samples": num_samples,
        "wells_std": wells_std,
        "wells_samp": wells_samp,
        "wells_neg": wells_neg,
        "wells_total": wells_total,
        "rb_samples": rb_samples,
        "rb_neg_ctrl": rb_neg_ctrl,
        "rb_standard": rb_standard,
        "rb_std_tube1": rb_std_tube1,
        "rb_std_dilutions": rb_std_dilutions,
        "rb_std_zero": rb_std_zero,
        "rb_dq_working": rb_dq_working,
        "rb_1x_subtotal": rb_1x_subtotal,
        "rb_1x_total": rb_1x_total,
        "rb_10x_needed": rb_10x_needed,
        "h2o_needed": h2o_needed,
        "dq_working_total": dq_working_total,
        "dq_stock_needed": dq_stock_needed,
        "dq_working_needed": dq_working_needed,
        "elastase_stock_needed": ELASTASE_STOCK_VOL,
        "lb_miller_needed": NEG_CTRL_LB_VOL,
        "sample_total_vol": SAMPLE_VOL_PER * num_samples,
        "concentrations": concentrations,
    }
