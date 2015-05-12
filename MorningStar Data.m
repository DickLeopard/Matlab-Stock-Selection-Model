function [varargout] = exlDatafromMorningstar(path_and_file,exlSheetName,exlRange)

% Access data through Excel
exlServerObj = actxserver('Excel.Application');
try 
    exlWkbkObj = exlServerObj.Workbooks;
    exlFile = exlWkbkObj.Open(path_and_file);
    exlSheet = exlFile.Sheets.Item(exlSheetName);
    exlRngObj = exlSheet.Range(exlRange);
    exlData = exlRngObj.Value;
catch
    disp(exception.message);  % what's the usage of this code? 
    exlServerObj.Quit;
    return
end
exlServerObj.Quit;
clear exlServerObj

% Find and remove all rows containing empty or invalid values in any column
% and create a cleaned set of data
notEmptyCells_Logical = ~cellfun('isempty',exlData);
notEmptyCells_Numerical = double(notEmptyCells_Logical);
rowsWithoutAnEmptyCell = prod(notEmptyCells_Numerical,2);
rowsIndex = find(rowsWithoutAnEmptyCell)';
trimmedExlData = cell([numel(rowsIndex),size(exlData,2)]);
count = 0;
for rowLoop = rowsIndex
    cellsWithoutText = ~cellfun('isclass',exlData(rowLoop,:),'char');
    cellsWithNumbersIndex = find(cellsWithoutText);
    invalidNumbers = isnan([exlData{rowLoop,cellsWithNumbersIndex}]);
    if ~any(invalidNumbers)
        count = count + 1;
        for columnLoop = 1:size(exlData,2)
            trimmedExlData{count,columnLoop} = exlData{rowLoop,columnLoop};
        end
    end
end
trimmedExlData = trimmedExlData([1:count],:);
clear rowsIndex rowLoop columnLoop count exlData

% Extract the cleaned set of output arguments
for kk = 1:nargout
    if ismember(kk,cellsWithNumbersIndex)
        varargout(kk) = {[trimmedExlData{:,kk}]'};
    else
        varargout(kk) = {trimmedExlData(:,kk)'};
    end
end
clear kk cellsWithNumbersIndex
