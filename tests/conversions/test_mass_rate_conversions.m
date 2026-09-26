% test that mass2rate and rate2mass are inverses of each other

delta_t = 0.1;
t = 2016:delta_t:2017;

% test from rate to mass back to rate
rate = rand(1,length(t));

[tmass,mass_est]=rate2mass(t,rate);
[trate,rate_est]=mass2rate(tmass,mass_est);

disp(['Rate: summed time error is ' num2str(sum(trate-t))]);
disp(['Rate: summed rate error is ' num2str(sum(rate_est-rate))]);

% test from mass to rate back to mass
mass = linspace(0,59,length(t));
[trate,rate_est]=mass2rate(t,mass);
[tmass,mass_est]=rate2mass(trate,rate_est);

disp(['Mass: summed time error is ' num2str(sum(tmass-t))]);
disp(['Mass: summed mass error is ' num2str(sum(mass_est-mass))]);