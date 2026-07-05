================================================================================
 PANEL STOCHASTIC FRONTIER (PSF) -- APPLICATION REPLICATION PACKAGE
================================================================================

This folder reproduces the empirical application results reported in the paper
by Tomioka, Yang and Zhang (2026), fourthcoming in the Journal of Business and 
Economics.
--------------------------------------------------------------------------------
 0. QUICK START -- HOW TO REPLICATE
--------------------------------------------------------------------------------

  1. Open MATLAB and set the *current/working folder* to this folder
     (the one containing this README.txt and main_application.m).

  2. Run the single driver script:

         >> main_application

     This is the ONLY file you run by hand. It adds the Functions/ and
     Scripts/ folders to the path, loads the data, executes the full
     estimation algorithm, computes standard errors, saves all output
     .mat files into Results/, and generates the figures.

  3. Outputs:
       - All numerical results are written to the Results/ folder as .mat files.
       - The figures are written to this folder as PDFs:
             Grouped_Frontiers_Application_SE.pdf   (frontiers with SE bands)
             Grouped_Frontiers_pihat.pdf            (estimated coefficient curves)
             ES_Application.pdf                     (efficiency / c-u summary)
       - Point estimates and standard errors (sigv, omegabar) are printed to
         the MATLAB Command Window.

--------------------------------------------------------------------------------
 1. FOLDER LAYOUT
--------------------------------------------------------------------------------

  Application/
  |
  |-- main_application.m      <-- RUN THIS. The driver for the whole pipeline.
  |-- Data_Application.xlsx       Panel data (N = 466 banks, T = 80 periods).
  |                               The code reads "Sheet1".
  |-- README.txt                  This file.
  |
  |-- Functions/                  Reusable functions called by the driver/scripts.
  |-- Scripts/                    Procedural scripts run in sequence by the driver.
  |-- Results/                    Destination for all saved .mat output files.
  |
  |-- SE_calculation.pdf          Derivation of the standard-error formulas.
  |-- PSF_TYZ_2026.pdf            Working paper version of the paper.
  |-- *.pdf (other)               Figures produced by the run (see above).

--------------------------------------------------------------------------------
 2. EXECUTION FLOW (what main_application.m does, in order)
--------------------------------------------------------------------------------

  main_application.m:
    - addpath('Functions'), addpath('Scripts'), addpath('Results')
    - Loads Data_Application.xlsx (Sheet5), takes logs, demeans/standardizes.
    - Sets the Params structure (T, N, m, p=5, K_group, initial values, bounds).
    - Calls the estimation engine:
         Est = Estimation(Params, y, x1..x5, x1_cl..x5_cl);
      Estimation.m (in Functions/) runs the whole algorithm -- SetStorage,
      Step 1 OLS sieve estimation, Steps 2-3 HAC grouping + information-criterion
      model selection (store_IC*), and Step 4/4' ML estimation of c/sigma_u and
      the mixture weight tau (loglik_K* / loglik_K*_prime via fminsearchbnd).
      It returns ALL variables it creates in the struct Est; the driver then
      unpacks Est into the workspace so the scripts below find what they expect.
    - Then runs, IN THIS ORDER (at the end of main_application.m):
         run getFrontiers              -> recovers fitted frontier curves;
                                          saves hat_* to Results/; internally
                                          calls SaveOutputs_Application to save
                                          and move the store_* .mat files.
         run Asymptotics               -> frontier SE covariance matrices (S_K*).
         run Frontiers_SE              -> assembles upper/lower SE bands.
         run Figures_Application       -> draws and prints the PDF figures.
         run sigv_omegabar_estimates   -> prints sigv and omegabar estimates + SEs.

--------------------------------------------------------------------------------
 3. FILE-BY-FILE DESCRIPTION
--------------------------------------------------------------------------------

ROOT
  main_application.m
      Master driver. Loads data, sets Params, runs the estimation algorithm,
      and calls the Scripts in sequence. Run this to replicate everything.
  Data_Application.xlsx
      Source data. Sheet1 is used (cost, x1, x2, y1, y2, y3); reshaped to a
      80 x 466 panel.

Functions/   (called as func(args); each returns values)
  Estimation.m          THE ESTIMATION ENGINE. Runs the full algorithm (Step 1
                        OLS sieve estimation, Step 2-3 HAC grouping + IC model
                        selection, Step 4/4' ML estimation of c/sigma_u/tau).
                        Called once by main_application.m; returns every variable
                        it creates in a struct, which the driver unpacks into the
                        workspace for the downstream scripts. (This is the code
                        that previously lived inline in main_application.m.)
  GenBasis.m            Builds the cosine sieve basis functions B0..B11.
  GenBm.m               Assembles the sieve matrix Bm (= B^m) for the intercept/
                        alpha terms given the degree m.
  GenxBm.m              Builds xBm = x_it (kron) basis, the regressor-sieve block.
  inputs_g.m            Constructs the per-group design array ZimLbar for group g
                        (used in the grouped likelihood / frontier steps).
  inputs_g0.m           Same as inputs_g.m but for the pooled (no-group) case.
  loglik_K1.m .. loglik_K4.m
                        NEGATIVE log-likelihood for the "no c-u group structure"
                        (single inefficiency component) case, for K = 1..4
                        sieve groups. Minimized in step 4.
  loglik_K1_prime.m .. loglik_K4_prime.m
                        NEGATIVE log-likelihood for the c-u MIXTURE case (the
                        "prime" / step-4' variant), K = 1..4. Parameter order is
                        omegabar = [c1, sqrt(sigu1), c2, sqrt(sigu2), tau].
  se_omegabar_K2_prime.m
                        Standard errors of the omegabar mixture parameters
                        (K=2, mixture case) by numerical (finite-difference)
                        observed information: SE = 1./sqrt(N * II). Output is
                        reordered to match the paper's display
                        [tau, c1, sigu1, c2, sigu2]; see the in-file comments and
                        wiki note application-02-omegabar-se-ordering.
  getFrontiers_K1_p5.m .. getFrontiers_K4_p5.m
                        Map estimated sieve coefficients (pi) back into the
                        fitted frontier curves alpha_hat and beta1..beta5_hat,
                        for p = 5 regressors and K = 1..4 groups.
  fminsearchbnd.m       Third-party bound-constrained Nelder-Mead optimizer
                        (wrapper around fminsearch). Used to maximize the
                        likelihoods subject to Params.LB/UB.

Scripts/   (run as scripts; they read/write variables in the base workspace)
  SetStorage.m              Preallocates the store_* cell/array containers.
  getFrontiers.m            Selects the chosen-K* pi estimates, calls the
                            getFrontiers_K*_p5 functions, saves hat_* curves to
                            Results/, and runs SaveOutputs_Application.
  Asymptotics.m             Builds the frontier SE covariance matrices (S_K*).
                            NOTE: the "SE for sigma_v" section here is only a
                            pointer comment -- the sigma_v SEs are computed inline
                            in main_application.m at each sigv2hat estimate.
  Frontiers_SE.m            Turns the S_K* matrices into pointwise upper/lower
                            standard-error bands for each frontier curve.
  Figures_Application.m     Loads the saved results and prints the PDF figures
                            (Grouped_Frontiers_Application_SE.pdf,
                            Grouped_Frontiers_pihat.pdf, ES_Application.pdf).
  sigv_omegabar_estimates.m Prints the sigma_v and omegabar estimates and their
                            standard errors, remapped to the paper's labeling.
  SaveOutputs_Application.m Saves the store_*/Params .mat files and moves them
                            into Results/ (via save + movefile). Called from
                            getFrontiers.m, not directly from the driver.

Results/   (output; created/overwritten by a run)
  Params.mat                  The Params structure used for the run.
  store_group_sieves.mat      Final HAC group assignment of firms.
  store_IC.mat, store_IC1_prime.mat, store_IC2_prime.mat
                              Information-criterion values used for model/mixture
                              selection.
  store_Kstar_L1.mat          Selected number of groups K* (benchmark c_lambda=1).
  store_MixD_L1.mat, store_MixD_L32.mat, store_MixD_L34.mat
                              Mixture-structure detection at c_tilde = 1, 3/2, 3/4.
  store_omegabar_K2.mat, store_sigvhat_K2.mat, store_pi_K2_g1/g2.mat
                              The selected-model parameter and coefficient estimates.
  hat_alphag*_Application.mat, hat_beta{1..5}g*_Application.mat
                              The fitted frontier curves per group, loaded by
                              Figures_Application.m.

--------------------------------------------------------------------------------
 4. CAUTIONS WHEN RUNNING THE MATLAB FILES
--------------------------------------------------------------------------------

  * Results/ MUST EXIST. SaveOutputs_Application.m moves files into Results/. If
    that folder is missing, movefile errors. It is included in this package; do
    not delete it. (Re-create an empty Results/ folder if you ever remove it.)

  * RE-RUNNING OVERWRITES. Each run overwrites the .mat files in Results/ and the
    output PDFs in the root. Back up any results you want to keep before re-running.
    If a target .mat is open in another program or locked, movefile/save will error.

  * RUN AS A WHOLE, IN ORDER. The Scripts are NOT standalone -- they read
    variables left in the base workspace by earlier steps and must run in the
    sequence the driver specifies (getFrontiers -> Asymptotics -> Frontiers_SE
    -> Figures_Application -> sigv_omegabar_estimates). Running a script on its
    own will error on undefined variables.

  * p IS FIXED AT 5. Params.p = 5 is assumed throughout: the frontier recovery
    functions are the *_p5 variants and the Asymptotics MB blocks are built for
    p = 5. Do not change Params.p without updating getFrontiers_K*_p5.m and
    Asymptotics.m accordingly.

  * MIXTURE LABEL SWITCHING. The two inefficiency mixture components are only
    identified up to relabeling. The optimizer may return the components in the
    opposite order to the paper, with tau reported as 1 - tau and
    (c1,sigu1) <-> (c2,sigu2) swapped. This is handled deliberately by the
    display remap in sigv_omegabar_estimates.m; the standard errors are reordered
    to match (see SE_calculation.pdf and the wiki note
    application-02-omegabar-se-ordering). Read the printed estimates with that
    remap in mind.

  * OPTIMIZER SENSITIVITY. Estimation uses Nelder-Mead (fminsearchbnd) subject to
    Params.LB/UB, started from the initial values Params.c1/sigu1/c2/sigu2/tau.
    Results can depend on these starting values and bounds; change them with care.

  * TOOLBOXES / DATA. Requires the Statistics and Machine Learning Toolbox
    (normcdf) and readtable (the .xlsx reader). Keep Data_Application.xlsx in
    this folder with its "Sheet1" tab intact.

  * RUNTIME. The estimation uses nested loops over N = 466 firms and T = 80
    periods and repeated likelihood evaluations; a full run can take a while.
    Let it finish before inspecting Results/.

================================================================================
