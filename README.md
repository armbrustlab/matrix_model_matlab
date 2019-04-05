# matrix_model_matlab
Prepping SeaFlow data and using it in the Hunter-Cevera 2-population matrix growth model.  You must have Kristen Hunter-Cevera's [phyto-division-rate-model](https://github.com/khuntercevera/phyto-division-rate-model) downloaded and in your Matlab path.

Matlab files include:
* `matrix_growth_rate_distributions.m`
    * Grabs PAR from the SQLite database, smooths it, and samples it at 10-minute intervals to prep for the matrix model.
    * Grabs VCT data from the `SeaFlow-OPP` directory and bins the volumes and carbon quotas.
* `PAR_2_matrix_growth.m`
    * Grabs PAR from the ship's underway data for cruises when it seems the sfl file doesn't have it.
* `matrix_growth_rate_driver.m`
    * Based on KHC's `call_to_opt_mvco.m`
    * Imports the PAR and size distributions saved during `matrix_growth_rate_distributions.m`.
    * Runs the paramter optimization routine (in parallel) for each day and stores parameters and growth rates to `modelresults`.  This is not yet exported because I have not yet had satisfactory convergence.
