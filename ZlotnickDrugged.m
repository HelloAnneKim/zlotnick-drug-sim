%% Mimick Zlotnick paper 1994 on Virus Capsid Assembly
m = InitializeModelMassActionAmount('Zlotnick1999');

%% 11 bimolecular reactions between monomer and species
%% monomer + n-species <-> (n+1)-species
m = AddCompartment(m, 'Solution', 3, 1);

%% Generate Drugged Capsids
drugged_capsids1 = GenerateDrugStringsOfN(1);
drugged_capsids2 = GenerateDrugStringsOfN(2);
drugged_capsids3 = GenerateDrugStringsOfN(3);
drugged_capsids4 = GenerateDrugStringsOfN(4);
drugged_capsids5 = GenerateDrugStringsOfN(5);
drugged_capsids6 = GenerateDrugStringsOfN(6);
drugged_capsids7 = GenerateDrugStringsOfN(7);
drugged_capsids8 = GenerateDrugStringsOfN(8);
drugged_capsids9 = GenerateDrugStringsOfN(9);
drugged_capsids10 = GenerateDrugStringsOfN(10);
drugged_capsids11 = GenerateDrugStringsOfN(11);
drugged_capsids12 = GenerateDrugStringsOfN(12);
all_drugged_capsids = cat(2,drugged_capsids1, drugged_capsids2, drugged_capsids3, drugged_capsids4,  drugged_capsids5, drugged_capsids6, drugged_capsids7, drugged_capsids8, drugged_capsids9, drugged_capsids10, drugged_capsids11, drugged_capsids12);

%% The monomers species and their original concentrations
drugless_capsids = {'n1','n2','n3','n4','n5','n6','n7','n8','n9','n10','n11','n12'};
m = AddState(m, 'n1', 'Solution', 2.5*10.^-6); % Molar (ref 1999 Zlotnick)
m = AddStatesWithSameConcentration(m, drugless_capsids(1,2:12), 'Solution', 0);

%% Add States of Drugged Capsids
m = AddState(m, 'd1', 'Solution', 2.5*10.^-6); % Molar 
% m = AddState(m, 'd1', 'Solution', 0*10.^-6); % Molar 
m = AddStatesWithSameConcentration(m, all_drugged_capsids(1,2:78), 'Solution', 0);

%% Outputs
m = AddOutputs(m, drugless_capsids, NamesToExpressions(drugless_capsids));
m = AddOutput(m, 'n12Times12', {'n12$' 12});
m = AddOutput(m, 'all_Subunits', {'n1$' 1; 'n2$' 2; 'n3$' 3; 'n4$' 4; 'n5$' 5; 'n6$' 6; 'n7$' 7; 'n8$' 8; 'n9$' 9; 'n10$' 10; 'n11$' 11; 'n12$' 12});
m = AddOutputs(m, all_drugged_capsids, NamesToExpressions(all_drugged_capsids));
m = AddOutput(m, 'n11d1Times12', {'n11d1$' 12});
m = AddOutput(m, 'n10d2Times12', {'n10d2$' 12});
m = AddOutput(m, 'n9d3Times12', {'n9d3$' 12});
m = AddOutput(m, 'n8d4Times12', {'n8d4$' 12});
m = AddOutput(m, 'n7d5Times12', {'n7d5$' 12});
m = AddOutput(m, 'n6d6Times12', {'n6d6$' 12});
m = AddOutput(m, 'n5d7Times12', {'n5d7$' 12});
m = AddOutput(m, 'n4d8Times12', {'n4d8$' 12});
m = AddOutput(m, 'n3d9Times12', {'n3d9$' 12});
m = AddOutput(m, 'n2d10Times12', {'n2d10$' 12});
m = AddOutput(m, 'n1d11Times12', {'n1d11$' 12});
m = AddOutput(m, 'd12Times12', {'d12$' 12});
multiplied12mers = {'n11d1Times12','n10d2Times12','n9d3Times12', 'n8d4Times12', 'n7d5Times12', 'n6d6Times12', 'n5d7Times12', 'n4d8Times12', 'n3d9Times12', 'n2d10Times12','n1d11Times12','d12Times12'};
all_outputs = cat(2, drugless_capsids, {'n12Times12'}, {'all_Subunits'}, all_drugged_capsids, multiplied12mers);

%% Parameters
k_fnuc = 10.^2.5 % M^-1*s^-1 (ref 1999 Zlotnick) dimer and trimer formation
m = AddParameter(m, 'k_fnuc', k_fnuc);
k_felong = 10.^4.5 % M^-1*s^-1 (ref 1999 Zlotnick)
m = AddParameter(m, 'k_felong', k_felong);
k_aContact = 10.^4.5; % M^-1 (ref 1999 Zlotnick)

S_2 = 5./2; % ratio of formation to dissociation (ref 1994 Zlotnick)
m = AddParameter(m, 'S_2', S_2);   
j_2 = 1; % number of contacts (ref 1994 Zlotnick as term 'c')
m = AddParameter(m, 'j_2', j_2);
k_a2 = k_aContact.^j_2.*S_2; % M^-1 (ref 1999 Zlotnick)
m = AddParameter(m, 'k_a2', k_a2);
k_b2 = k_felong./k_a2 % s^-1 (ref 1999 Zlotnick)
m = AddParameter(m, 'k_b2', k_b2);

S_3 = 2./3;
m = AddParameter(m, 'S_3', S_3);
j_3 = 2;
m = AddParameter(m, 'j_3', j_3);
k_a3 = k_aContact.^j_3.*S_3;
m = AddParameter(m, 'k_a3', k_a3);
k_b3 = k_felong./k_a3
m = AddParameter(m, 'k_b3', k_b3);

S_4 = 3./2;
m = AddParameter(m, 'S_4', S_4);
j_4 = 2;
m = AddParameter(m, 'j_4', j_4);
k_a4 = k_aContact.^j_4.*S_4;
m = AddParameter(m, 'k_a4', k_a4);
k_b4 = k_felong./k_a4;
m = AddParameter(m, 'k_b4', k_b4);

S_5 = 4./2;
m = AddParameter(m, 'S_5', S_5);
j_5 = 2;
m = AddParameter(m, 'j_5', j_5);
k_a5 = k_aContact.^j_5.*S_5;
m = AddParameter(m, 'k_a5', k_a5);
k_b5 = k_felong./k_a5;
m = AddParameter(m, 'k_b5', k_b5);

S_6 = 1./5;
m = AddParameter(m, 'S_6', S_6);
j_6 = 3;
m = AddParameter(m, 'j_6', j_6);
k_a6 = k_aContact.^j_6.*S_6;
m = AddParameter(m, 'k_a6', k_a6);
k_b6 = k_felong./k_a6;
m = AddParameter(m, 'k_b6', k_b6);

S_7 = 5./1;
m = AddParameter(m, 'S_7', S_7);
j_7 = 2;
m = AddParameter(m, 'j_7', j_7);
k_a7 = k_aContact.^j_7.*S_7;
m = AddParameter(m, 'k_a7', k_a7);
k_b7 = k_felong./k_a7;
m = AddParameter(m, 'k_b7', k_b7);

S_8 = 2./4;
m = AddParameter(m, 'S_8', S_8);
j_8 = 3;
m = AddParameter(m, 'j_8', j_8);
k_a8 = k_aContact.^j_8.*S_8;
m = AddParameter(m, 'k_a8', k_a8);
k_b8 = k_felong./k_a8;
m = AddParameter(m, 'k_b8', k_b8);

S_9 = 2./3;
m = AddParameter(m, 'S_9', S_9);
j_9 = 3;
m = AddParameter(m, 'j_9', j_9);
k_a9 = k_aContact.^j_9.*S_9;
m = AddParameter(m, 'k_a9', k_a9);
k_b9 = k_felong./k_a9;
m = AddParameter(m, 'k_b9', k_b9);

S_10 = 3./2;
m = AddParameter(m, 'S_10', S_10);
j_10 = 3;
m = AddParameter(m, 'j_10', j_10);
k_a10 = k_aContact.^j_10.*S_10;
m = AddParameter(m, 'k_a10', k_a10);
k_b10 = k_felong./k_a10;
m = AddParameter(m, 'k_b10', k_b10);

S_11 = 2./5;
m = AddParameter(m, 'S_11', S_11);
j_11 = 4;
m = AddParameter(m, 'j_11', j_11);
k_a11 = k_aContact.^j_11.*S_11;
m = AddParameter(m, 'k_a11', k_a11);
k_b11 = k_felong./k_a11;
m = AddParameter(m, 'k_b11', k_b11);

S_12 = 1./12;
m = AddParameter(m, 'S_12', S_12);
j_12 = 5;
m = AddParameter(m, 'j_12', j_12);
k_a12 = k_aContact.^j_12.*S_12;
m = AddParameter(m, 'k_a12', k_a12);
k_b12 = k_felong./k_a12;
m = AddParameter(m, 'k_b12', k_b12);

%% 11 BiMolecular Reactions 
m = AddReaction(m, '', {'n1', 'n1'}, {'n2'}, 'k_fnuc', 'k_b2'); % n1+n1 -> n2 w/ k_f forward and k_b2 backwards
m = AddReaction(m, '', {'n1', 'n2'}, {'n3'}, 'k_felong', 'k_b3');
m = AddReaction(m, '', {'n1', 'n3'}, {'n4'}, 'k_felong', 'k_b4');
m = AddReaction(m, '', {'n1', 'n4'}, {'n5'}, 'k_felong', 'k_b5');
m = AddReaction(m, '', {'n1', 'n5'}, {'n6'}, 'k_felong', 'k_b6');
m = AddReaction(m, '', {'n1', 'n6'}, {'n7'}, 'k_felong', 'k_b7');
m = AddReaction(m, '', {'n1', 'n7'}, {'n8'}, 'k_felong', 'k_b8');
m = AddReaction(m, '', {'n1', 'n8'}, {'n9'}, 'k_felong', 'k_b9');
m = AddReaction(m, '', {'n1', 'n9'}, {'n10'}, 'k_felong', 'k_b10');
m = AddReaction(m, '', {'n1', 'n10'}, {'n11'}, 'k_felong', 'k_b11');
m = AddReaction(m, '', {'n1', 'n11'}, {'n12'}, 'k_felong', 'k_b12');

%% Drug Reactions
m = AddReaction(m, '', {'n1', 'd1'}, {'n1d1'}, 'k_fnuc', 'k_b2');
m = AddReaction(m, '', {'d1', 'd1'}, {'d2'}, 'k_fnuc', 'k_b2');
m = AddDrugReactionsN(m, '', 3, 'k_felong', 'k_b3');
m = AddDrugReactionsN(m, '', 4, 'k_felong', 'k_b4');
m = AddDrugReactionsN(m, '', 5, 'k_felong', 'k_b5');
m = AddDrugReactionsN(m, '', 6, 'k_felong', 'k_b6');
m = AddDrugReactionsN(m, '', 7, 'k_felong', 'k_b7');
m = AddDrugReactionsN(m, '', 8, 'k_felong', 'k_b8');
m = AddDrugReactionsN(m, '', 9, 'k_felong', 'k_b9');
m = AddDrugReactionsN(m, '', 10, 'k_felong', 'k_b10');
m = AddDrugReactionsN(m, '', 11, 'k_felong', 'k_b11');
m = AddDrugReactionsN(m, '', 12, 'k_felong', 'k_b12');

m = FinalizeModel(m);

%% Construct experiment
con = experimentInitialValue(m, [], [], [], 'InitialValueExperiment');

%% Simulate
tF = 600; % final time in seconds
sim1 = SimulateSystem(m, con, tF); %this is where all the data is

figure
indices_of_output = [1,12,13,15, 92:104];
times = linspace(0, tF, 100);
% plot(times, sim1.y(times))
plot(times, sim1.y(times, indices_of_output))
% legend(cat(2, drugless_capsids, {'n12Times12'}, 'n2d1'))
legend(all_outputs(indices_of_output));
xlabel('Time')
ylabel('Amount')
title('Drugging Zlotnick Model at 2.5*10.^-6 drug');


