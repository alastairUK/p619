classdef P619Tests < matlab.unittest.TestCase
    % P619Tests
    % Unit/regression tests for P619.
    properties
        ITURP619
    end

    methods (TestClassSetup)
        function createP619Instance(testCase)
            % createP619Instance
            % Construct the P619 instance once for the test class.
            testCase.ITURP619 = P619;
        end
    end

    methods (Test)
        function tl_free_space(testCase)
            % tl_free_space
            % Validate free-space loss against expected values..
            T = readtable(P619Tests.dataFile('tl_free_space_cases.csv'));
            absTol = 1e-6;

            for i = 1:height(T)
                actual = testCase.ITURP619.tl_free_space(T.fGHz(i), T.dkm(i));
                msg = sprintf('Test failed for fGHz=%.6g, dkm=%.6g', T.fGHz(i), T.dkm(i));
                testCase.verifyEqual(actual, T.expected(i), 'AbsTol', absTol, msg);
            end
        end

        function co_and_cross_polar_attenuation(testCase)
            % co_and_cross_polar_attenuation
            % Validate co- and cross-polar attenuation outputs.
            T = readtable(P619Tests.dataFile('co_and_cross_polar_attenuation.csv'));
            absTol = 1e-6;

            for k = 1:height(T)
                xpd = T.XPD(k);
                [Across, Aco] = testCase.ITURP619.co_and_cross_polar_attenuation(xpd);

                msg = sprintf('Mismatch at XPD = %.6g dB', xpd);
                testCase.verifyEqual(Aco,    T.Aco_expected(k),    'AbsTol', absTol, msg);
                testCase.verifyEqual(Across, T.Across_expected(k), 'AbsTol', absTol, msg);
            end
        end

        function beam_spreading_loss_834_8(testCase)
            % beam_spreading_loss_834_8
            % Regression test for -8 beam spreading behaviour (|10*log10(B)|).

            obj = testCase.ITURP619;
            T = readtable(P619Tests.dataFile('beam_spreading_cases.csv'));
            absTol = 1e-9;

            for i = 1:height(T)
                actVal = obj.beam_spreading_loss(T.theta0_deg(i), T.h_km(i));
                testCase.verifyEqual(actVal, T.expected_p834_8(i), 'AbsTol', absTol);
            end
        end

        function beam_spreading_loss_834_9(testCase)
            % beam_spreading_los_834_9
            % Validate newer  definition (Abs = -10*log10(B)).

            obj = testCase.ITURP619;
            T = readtable(P619Tests.dataFile('beam_spreading_cases.csv'));
            absTol = 1e-9;

            for i = 1:height(T)
                actVal = obj.beam_spreading_loss_834_9(T.theta0_deg(i), T.h_km(i));
                testCase.verifyEqual(actVal, T.expected_p834_9(i), 'AbsTol', absTol);
            end
        end

        function p676d11_ga_specific_attenuation_validation(testCase)
            % Test: p676d11_ga_specific_attenuation_validation
            % Validate P.676-11 specific attenuation components returned by p676d11_ga()
            % (oxygen g0 and water vapour gw) against ITU validation sheet values.

            csvPath = fullfile(fileparts(mfilename('fullpath')), "testdata", ...
                "p676d11_ga_validation_specific_attenuation.csv");

            tbl = readtable(csvPath);

            req = ["f_GHz","p_hPa","rho_gm3","T_K", ...
                "expected_gamma0_dB_per_km","expected_gammaw_dB_per_km"];
            testCase.verifyTrue(all(ismember(req, string(tbl.Properties.VariableNames))), ...
                "CSV missing required columns.");

            absTol = 1e-8;

            for k = 1:height(tbl)
                f = tbl.f_GHz(k);
                p = tbl.p_hPa(k);
                rho = tbl.rho_gm3(k);
                T = tbl.T_K(k);

                exp_g0 = tbl.expected_gamma0_dB_per_km(k);
                exp_gw = tbl.expected_gammaw_dB_per_km(k);

                [g0, gw] = testCase.ITURP619.p676d11_ga(f, p, rho, T);

                msg = sprintf("Row %d: f=%.6g GHz, pdry=%.6g hPa, rho=%.6g g/m3, T=%.6g K", ...
                    k, f, p, rho, T);

                testCase.verifyEqual(g0, exp_g0, "AbsTol", absTol, msg + " (g0)");
                testCase.verifyEqual(gw, exp_gw, "AbsTol", absTol, msg + " (gw)");
            end
        end

        function p618_xpd_scaling_validation(testCase)
            % Test: p618_xpd_scaling_validation
            % Validate scaling of long-term XPD between frequencies/polarizations.
            T = readtable(P619Tests.dataFile('p618_xpd_scaling_cases.csv'));

            absTol = 1e-9;
            for k = 1:height(T)
                act = testCase.ITURP619.p618_xpd_scaling(T.xpd1_dB(k), T.f1_GHz(k), T.tau1_deg(k), T.f2_GHz(k), T.tau2_deg(k));
                msg = sprintf("Row %d: xpd1=%.6g dB, f1=%.6g GHz, tau1=%.6g deg, f2=%.6g GHz, tau2=%.6g deg", ...
                    k, T.xpd1_dB(k), T.f1_GHz(k), T.tau1_deg(k), T.f2_GHz(k), T.tau2_deg(k));
                testCase.verifyEqual(act, T.expected_xpd2_dB(k), "AbsTol", absTol, msg);
            end
        end

        function p618_hydrometeor_xpd_validation(testCase)
            % Test: p618_hydrometeor_xpd_validation
            % Validate P.618-14 cross-polarization discrimination
            % from rain statistics
            csvPath = fullfile(fileparts(mfilename('fullpath')), "testdata", ...
                "P618_14_XPD.csv");

            tbl = readtable(csvPath);

            req = ["f_GHz","P_percent","A_p","ElevationAngle_deg","PolarizationTiltAngle_deg","XPD_p"];
            testCase.verifyTrue(all(ismember(req, string(tbl.Properties.VariableNames))), ...
                "CSV missing required columns.");

            absTol = 1e-8;

            for k = 1:height(tbl)
                f = tbl.f_GHz(k);
                p = tbl.P_percent(k);
                A_p = tbl.A_p(k);
                ElevationAngle_deg = tbl.ElevationAngle_deg(k);
                PolarizationTiltAngle_deg = tbl.PolarizationTiltAngle_deg(k);

                if (ElevationAngle_deg > 60)
                    warning("Row %d: ElevationAngle_deg=%.6g exceeds model validity range (<=60 deg). Skipping.", ...
                        k, ElevationAngle_deg);
                    continue;
                end
                
                exp_XPD_p = tbl.XPD_p(k);

                [xpd] = testCase.ITURP619.p618_hydrometeor_xpd(f, p, A_p, PolarizationTiltAngle_deg, ElevationAngle_deg);
                msg = sprintf("Row %d: f=%.6g GHz, P=%.6g percent, A_p=%.6g, ElevationAngle_deg=%.6g, PolarizationTiltAngle_deg=%.6g", ...
                    k, f, p, A_p, ElevationAngle_deg, PolarizationTiltAngle_deg);

                testCase.verifyEqual(xpd, exp_XPD_p, "AbsTol", absTol, msg + " (XPD_p)");
            end
        end
        function beam_spreading_loss_834_9_random_approval(testCase)
            % beam_spreading_loss_834_9_random_approval
            % Approval-style CSV test for beam_spreading_loss_834_9().
            % - Generates deterministic random inputs (theta0_deg, h_km)
            % - Computes p834_9 results
            % - If approved CSV exists: compares against it
            % - If not: creates approved CSV
            % Files written under ./test_data next to this test file:
            %   beam_spreading_random_p834_9_approved.csv
            %   beam_spreading_random_p834_9_current.csv
        
            obj = testCase.ITURP619;
        
            seed = 12345;
            N = 1000;
        
            rng(seed, "twister");
            theta0_deg = 0.001 + (9.999 - 0.001) * rand(N, 1);
            h_km       = 0.0   + (4.999 - 0.0)   * rand(N, 1);
        
            p834_9 = zeros(N, 1);
            for i = 1:N
                p834_9(i) = obj.beam_spreading_loss_834_9(theta0_deg(i), h_km(i));
            end
        
            testDir = fileparts(which("P619Tests"));
            outDir = fullfile(testDir, "testdata");
            if exist(outDir, "dir") ~= 7
                mkdir(outDir);
            end
        
            approvedFile = fullfile(outDir, "beam_spreading_random_p834_9_approved.csv");
            currentFile  = fullfile(outDir, "beam_spreading_random_p834_9_current.csv");
        
            Tcur = table(theta0_deg, h_km, p834_9);
            writetable(Tcur, currentFile);
        
            if exist(approvedFile, "file") ~= 2
                writetable(Tcur, approvedFile);
                testCase.verifyTrue(true);
                return;
            end
        
            Tapp = readtable(approvedFile);
        
            req = ["theta0_deg","h_km","p834_9"];
            testCase.verifyTrue(all(ismember(req, string(Tapp.Properties.VariableNames))), ...
                "Approved CSV missing required columns.");
        
            testCase.verifyEqual(height(Tapp), height(Tcur), "Approved/current row count mismatch.");
        
            absTol = 1e-12;
            testCase.verifyEqual(Tcur.theta0_deg, Tapp.theta0_deg, "AbsTol", absTol, "theta0_deg mismatch.");
            testCase.verifyEqual(Tcur.h_km,       Tapp.h_km,       "AbsTol", absTol, "h_km mismatch.");
            testCase.verifyEqual(Tcur.p834_9,     Tapp.p834_9,     "AbsTol", absTol, "p834_9 mismatch.");
        end

        function get_interp2_Nwet_Annual_time_location_validation(testCase)
            % Test: get_interp2_Nwet_Annual_time_location_validation
            % Validates get_interp2_Nwet_Annual_time_location() against known values.
            csvPath = fullfile(fileparts(mfilename('fullpath')), "testdata", ...
                "P618_14_A_Scint.csv");

            tbl = readtable(csvPath);

            req = ["Lat_degN","Lon_degE","MedianNwet_P453"];
            testCase.verifyTrue(all(ismember(req, string(tbl.Properties.VariableNames))), ...
                "CSV missing required columns.");

            absTol = 1e-7;

            for k = 1:height(tbl)
                Lat_degN = tbl.Lat_degN(k);
                Lon_degE = tbl.Lon_degE(k);

                exp_MedianNwet = tbl.MedianNwet_P453(k);

                act_MedianNwet = testCase.ITURP619.get_interp2_Nwet_Annual_time_location(50.0, Lon_degE, Lat_degN);

                msg = sprintf("Row %d: Lat_degN=%.6g, Lon_degE=%.6g", ...
                    k, Lat_degN, Lon_degE);

                testCase.verifyEqual(act_MedianNwet, exp_MedianNwet, "AbsTol", absTol, msg);
            end
        end

        function tropospheric_scintillation_validation(testCase)
            % Test: tropospheric_scintillation_validation
            % Validates tropospheric_scintillation() against known values.
            csvPath = fullfile(fileparts(mfilename('fullpath')), "testdata", ...
                "P618_14_A_Scint_with_P619_Variation.csv");

            tbl = readtable(csvPath);

            req = ["f_GHz","ElevationAngle_deg","AntennaDiameter_m", ...
                   "AntennaEfficiency_h","p_percent","MedianNwet_P453","Ast_dB"];
            testCase.verifyTrue(all(ismember(req, string(tbl.Properties.VariableNames))), ...
                "CSV missing required columns.");

            absTol = 1e-7;

            for k = 1:height(tbl)
                f = tbl.f_GHz(k);
                ElevationAngle_deg = tbl.ElevationAngle_deg(k);
                AntennaDiameter_m = tbl.AntennaDiameter_m(k);
                AntennaEfficiency_h = tbl.AntennaEfficiency_h(k);
                p_percent = tbl.p_percent(k);
                MedianNwet_P453 = tbl.MedianNwet_P453(k);

                exp_Ast_dB = tbl.Ast_dB(k);

                Ga = 10.0*log10(AntennaEfficiency_h*((pi*AntennaDiameter_m*f*1e9)/3e8)^2);
                act_Ast_dB = testCase.ITURP619.tropospheric_scintillation(f, p_percent, MedianNwet_P453, ElevationAngle_deg, Ga);

                msg = sprintf("Row %d: f=%.6g, ElevationAngle_deg=%.6g", ...
                    k, f, ElevationAngle_deg);

                testCase.verifyEqual(act_Ast_dB, exp_Ast_dB, "AbsTol", absTol, msg);
            end
        end

        function p452_path_profile_approval(testCase)
            % p452_path_profile_approval
            % Generates deterministic random path profiles and calculates
            % horizon parameters using P619.p452_path_profile.
            %
            % Outputs a CSV compatible with the C++ Rec452Tests reader:
            % ID, h_ts, ae, Expected_Theta, Expected_Dist, NumPoints, d0, h0, d1, h1...
            
            obj = testCase.ITURP619;
            
            % 1. Setup Deterministic Randomness
            seed = 54321;
            rng(seed, "twister");
            
            num_tests = 500;
            
            % 2. Setup File Paths
            testDir = fileparts(which('P619Tests'));
            outDir = fullfile(testDir, 'testdata');
            if exist(outDir, 'dir') ~= 7
                mkdir(outDir);
            end
            
            approvedFile = fullfile(outDir, 'p452_horizon_approved.csv');
            currentFile  = fullfile(outDir, 'p452_horizon_current.csv');
            
            % 3. Generate Data and Write to Current File
            fid = fopen(currentFile, 'w');
            % Use onCleanup to ensure file closes even if test fails midway
            c = onCleanup(@() fclose(fid));
            
            % Header matching C++ reader expectation
            fprintf(fid, 'ID,h_ts,ae,Expected_Theta,Expected_Dist,NumPoints,ProfileData\n');
            
            for i = 1:num_tests
                % --- Generate Inputs ---
                
                % Number of points (between 5 and 100)
                n = randi([5, 100]);
                
                % Total distance (10km to 200km)
                d_total = 10 + (190 * rand());
                
                % Create distance vector (uniform spacing)
                d = linspace(0, d_total, n);
                
                % Create height vector (random terrain 0-500m with a random hill)
                % This ensures the horizon isn't always the last point
                base_h = 500 * rand(1, n);
                hill_center = randi([2, n-1]);
                hill_height = 200 * rand();
                % Gaussian hill
                h = base_h + hill_height * exp(-((1:n) - hill_center).^2 / 10);
                
                % Tx Height (amsl)
                h_ts = h(1) + 10 + (50 * rand()); 
                
                % Effective Earth Radius (standard to super-refractive)
                ae = 6371 * (1 + 3 * rand()); 
                
                % --- Calculate (The Oracle) ---
                [d_hoz, th_hoz] = obj.p452_path_profile(d, h, h_ts, ae);
                
                % Handle potential empty returns (though unlikely with this generation)
                if isempty(d_hoz), d_hoz = 0; end
                if isempty(th_hoz), th_hoz = 0; end
                
                % --- Write Row ---
                % Fixed columns
                fprintf(fid, '%d,%.6f,%.6f,%.6f,%.6f,%d', i, h_ts, ae, th_hoz, d_hoz, n);
                
                % Variable profile columns (d, h pairs)
                for k = 1:n
                    fprintf(fid, ',%.6f,%.6f', d(k), h(k));
                end
                fprintf(fid, '\n');
            end
            
            % Force close before reading back
            clear c; 
            
            % 4. Approval Logic
            if exist(approvedFile, 'file') ~= 2
                % If approved file doesn't exist, create it from current
                copyfile(currentFile, approvedFile);
                warning('Approved file did not exist. Created: %s', approvedFile);
                testCase.verifyTrue(true);
                return;
            end
            
            % Read both files as text to ensure exact match
            currentText = fileread(currentFile);
            approvedText = fileread(approvedFile);
            
            % Verify
            testCase.verifyEqual(currentText, approvedText, ...
                'Generated horizon data does not match the approved file. Logic may have changed.');
        end


        function p835_reference_atmosphere_validation(testCase)
            testCase.p835_profile_validation_helper('p835_reference_atmosphere.csv', 1);
        end
        function p835_low_latitude_annual_reference_validation(testCase)
            testCase.p835_profile_validation_helper('p835_low_latitude_annual_reference.csv', 2);
        end
        function p835_mid_latitude_summer_reference_validation(testCase)
            testCase.p835_profile_validation_helper('p835_mid_latitude_summer_reference.csv', 31);
        end
        function p835_mid_latitude_winter_reference_validation(testCase)
            testCase.p835_profile_validation_helper('p835_mid_latitude_winter_reference.csv', 32);
        end
        function p835_high_latitude_summer_reference_validation(testCase)
            testCase.p835_profile_validation_helper('p835_high_latitude_summer_reference.csv', 41);
        end
        function p835_high_latitude_winter_reference_validation(testCase)
            testCase.p835_profile_validation_helper('p835_high_latitude_winter_reference.csv', 42);
        end
    end

    methods (Access = private)
        function p835_profile_validation_helper(testCase, csvFile, profileType)
            % Helper to validate p835 profile against known values in a CSV.
            csvPath = fullfile(fileparts(mfilename('fullpath')), "testdata", csvFile);
            tbl = readtable(csvPath);
            req = ["height_km","temp_K","pressure_hPa","wvd_gm3","n"];
            testCase.verifyTrue(all(ismember(req, string(tbl.Properties.VariableNames))), ...
                "CSV missing required columns.");
            absTol = 1e-7;
            for k = 1:height(tbl)
                height_km = tbl.height_km(k);
                exp_temp_K = tbl.temp_K(k);
                exp_pressure_hPa = tbl.pressure_hPa(k);
                exp_wvd_gm3 = tbl.wvd_gm3(k);
                exp_n = tbl.n(k);
                [act_temp_K, act_pressure_hPa, act_wvd_gm3, act_n] = ...
                    testCase.ITURP619.p835_std_atm_profiles(height_km, profileType);
                testCase.verifyEqual(act_temp_K,      exp_temp_K,      "AbsTol", absTol, ...
                    sprintf("Row %d: height_km=%.6g (temp_K)", k, height_km));
                testCase.verifyEqual(act_pressure_hPa, exp_pressure_hPa, "AbsTol", absTol, ...
                    sprintf("Row %d: height_km=%.6g (pressure_hPa)",k, height_km));
                testCase.verifyEqual(act_wvd_gm3,     exp_wvd_gm3,     "AbsTol", absTol, ...
                    sprintf("Row %d: height_km=%.6g (wvd_gm3)", k, height_km));
                testCase.verifyEqual(act_n,           exp_n,           "AbsTol", absTol, ...
                    sprintf("Row %d: height_km=%.6g (n)", k, height_km));
            end
        end
    end

    methods (Static, Access = private)
        function p = dataFile(name)
            % dataFile
            % Build an absolute path to a CSV under ./testdata next to this file.
            % Inputs:
            %   name - string/char, CSV filename (e.g. 'my_test.csv')
            % Outputs:
            %   p - absolute path to the CSV
            testDir = fileparts(which('P619Tests'));
            p = fullfile(testDir, 'testdata', name);
        end
    end
end
