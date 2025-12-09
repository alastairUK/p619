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
