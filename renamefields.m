function new_struct = renamefields(old_struct, oldnames, newnames, ~) % ~ noerror
% 2017-07-09  Matlab2006+  Copyright (c) 2017, W J Whiten  BSD License

if nargin == 4
    noerror = true;
else
    noerror = false;
end
% check and adjust inputs
if ~isstruct(old_struct)
    error('changefields:  First argument must be a struct')
end
if ischar(oldnames)
    oldnames = {oldnames};
end
if ischar(newnames)
    newnames = {newnames};
end
if length(oldnames) ~= length(newnames)
    error('changefields:  Number of names not equal');
end

% undo struct
names = fieldnames(old_struct);
names1 = names;
values = struct2cell(old_struct);

% change names
for i = 1:length(oldnames)
    ind = find(strcmp(oldnames{i}, names));
    if isempty(ind) && ~noerror
        error(['changefields:  Name  ''', oldnames{i}, '''  not in struct']);
    elseif ~isempty(ind)
        names1{ind} = newnames{i};
    end
end

% create new struct
new_struct = cell2struct(values, names1);

end