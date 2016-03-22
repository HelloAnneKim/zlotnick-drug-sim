function model = AddOutputs( model, names, expressions)
%AddOutputs Add multiple outputs
%   AddOutput()
[m,n] = size(names);
for i = 1:n
    name = char(names(1,i));
    expression = expressions(1,i);
    model = AddOutput(model, name, expression);
end
end

