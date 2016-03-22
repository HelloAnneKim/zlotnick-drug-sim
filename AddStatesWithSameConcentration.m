function m = AddStatesWithSameConcentration(m, names, compartment, seed)
%AddStatesWithSameConcentration Add multiple state species using AddState
%   AddState()

for name = names
    m = AddState(m, char(name{1}), compartment, seed);
end
end

