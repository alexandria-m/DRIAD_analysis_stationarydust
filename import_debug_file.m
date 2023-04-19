var_set = 0;
sigma_set = 0;
fid = fopen([path folder dataset{d} name '_debug.txt']);
if(fid ~= -1)
   fprintf(['\n i = ' num2str(d)])
   fprintf(['; importing ' dataset{d} '\n']);
   while ~feof(fid)
       tline = fgetl(fid);
       if(var_set == 1)
          [~,q] = find(tline == ' ');
          tempvar = deblank(tline(1:q(1)-1));
          tempval = tline(q(end)+1:end);
          eval([tempvar ' = ' tempval ';']);
            % stop reading and saving variables when you get to these
            if(tempvar == "BOLTZMANN" || ...
            	tempvar == "TIME_EVOL" || ...
            	tempvar == "E_FIELDR" || ...
            	tempvar == "NUM_DEN_GAS" || ...
            	tempvar == "dx")
            	var_set = 0;
             end
        end
        if(sigma_set > 0)
          switch sigma_set
           case 1
              sigma_set = sigma_set + 1;
           case 2
              first_10_sigma_i_tot = str2double(regexp(tline,' ','split'));
              first_10_sigma_i_tot(end) = [];
              sigma_set = sigma_set + 1;
           case 3
              sigma_set = sigma_set + 1;
           case 4
              first_10_sigma_i1 = str2double(regexp(tline,' ','split'));
              first_10_sigma_i1(end) = [];
              sigma_set = sigma_set + 1;
           case 5
              sigma_set = sigma_set + 1;
           case 6
              first_10_sigma_i2 = str2double(regexp(tline,' ','split'));
              first_10_sigma_i2(end) = [];
              sigma_set = sigma_set + 1;
           case 7
              sigma_set = sigma_set + 1;
           case 8
              [~,q] = find(tline == ' ');
              tempvar = deblank(tline(1:q(1)-2));
              tempval = tline(q(end)+1:end);
              eval([tempvar ' = ' tempval ';']);
              sigma_set = sigma_set + 1;
          end
        end
        % start reading and saving the variables after these lines
        if(tline == "-- Constants --" || ...
            tline == "-- User Parameters --" || ...
            tline == "-- Derived Parameters --" || ...
            tline == "-- Super Ion Parameters --" || ...
            tline == "-- Further Derived Parameters --")
            var_set = 1;
        end
        if(tline == "-- Sigma Values --")
            sigma_set = 1;
        end
        if(tline == "-- Time Step Commands --")
            sigma_set = 0;
        end
    end
end