function m = AddDrugReactionsN(m, name, N, kForward, kReverse)
%AddDrugReactionsN adds the drug reactions for a given N number of subunits
%in each drugged molecule
%   m = AddDrugReactionsN(m, name, N, kForward, kReverse)
% 
%SO IT ASSUMES THAT ALL THE DRUG REACTIONS HAVE THE SAME NAME AND KFORWARDS
%AND REVERSES
%
%For example, m = AddReaction(m, '', 3, 'k_f', 'k_b2');
% adds reactions:
% ['n2d1', 'n1d2', 'd3']
% m = AddReaction(m, '', {'d1', 'n2'}, {'n2d1'}, 'k_f', 'k_b2');
% m = AddReaction(m, '', {'n1', 'n1d1'}, {'n2d1'}, 'k_f', 'k_b2');

% m = AddReaction(m, '', {'d1', 'n1d1'}, {'n1d2'}, 'k_f', 'k_b2');
% m = AddReaction(m, '', {'n1', 'd2'}, {'n1d2'}, 'k_f', 'k_b2');

% m = AddReaction(m, '', {'d1', 'd2'}, {'d3'}, 'k_f', 'k_b2');
n1_reactions = 0;
d1_reactions = 0;
for drug = 1:N
    n = N - drug; % amount of non-drug subunits
    % combo subunit and drug
    if drug<N
        drugged = strcat('n',num2str(n),'d',num2str(drug));
        % drug monomer reaction
        if drug==1
           without_d1 = strcat('n',num2str(n));
           d1_reactions = d1_reactions+1; 
           m = AddReaction(m, name, {'d1', without_d1}, {drugged}, kForward, kReverse);
        else
           without_d1 = strcat('n',num2str(n),'d',num2str(drug-1));
           d1_reactions = d1_reactions+1; 
           m = AddReaction(m, name, {'d1', without_d1}, {drugged}, kForward, kReverse);
        end
        % subunit monomer reaction
        if n==1
            without_n1 = strcat('d',num2str(drug));
            n1_reactions = n1_reactions+1;
            m = AddReaction(m, name, {'n1', without_n1}, {drugged}, kForward, kReverse);
        else
            without_n1 = strcat('n',num2str(n-1),'d',num2str(drug));
            n1_reactions = n1_reactions+1;
            m = AddReaction(m, name, {'n1', without_n1}, {drugged}, kForward, kReverse);
        end
    % just drug monomers and subunits
    else
        drugged = strcat('d',num2str(drug));
        without_d1 = strcat('d',num2str(drug-1));
        d1_reactions = d1_reactions+1; 
        m = AddReaction(m, name, {'d1', without_d1}, {drugged}, kForward, kReverse);
    end
end
end

