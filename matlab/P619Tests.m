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
