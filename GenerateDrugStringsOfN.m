function list_of_drug_strings_for_N = GenerateDrugStringsOfN(N)
%GenerateDrugStringsOfN generates a list of drug strings
%
% list_of_drug_strings_for_N = GenerateDrugStringsOfN(N)
%
% For example, GenerateDrugStringsOfN(1) = {'d1'}; GenerateDrugStringsOfN(3) = {'n2n1', 'n1d2', 'd3'}
drugged_molecules = {};
for drug = 1:N
    if drug<N
        drugged = {strcat('n',num2str(N-drug),'d',num2str(drug))};
    else
        drugged = {strcat('d',num2str(drug))};
    end
    drugged_molecules = cat(2, drugged_molecules, drugged);
end
list_of_drug_strings_for_N = drugged_molecules;
end