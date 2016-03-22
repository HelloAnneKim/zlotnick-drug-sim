function expressions = NamesToExpressions( names )
%NamesToExpressions turns an array of names to an array of expressions
%   list_of_expressions = NamesToExpressions( list_of_names )
% For example, {'n1$','n2$'} = NamesToExpressions( {'n1', 'n2'} )
[m,n] = size(names);
expressions = {};
for i = 1:n
    name = char(names(1,i));
    expression = {strcat(name, '$')};
    expressions = cat(2, expressions, expression);
end
end

