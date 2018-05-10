function create_heatmap(cfmat_percent, cfmat_rawnum, fieldnames_beh, indicator,figure_indicator)

% CASES
% 1. cross validation (same training and testing set) for percents
% 2. case 1 but with raw numbers
% 3. test rat is not in training set for percents
% 4. case 3 but with raw numbers

switch indicator 
    case 1 
        cfmat = cfmat_percent ; 
        indices = all(~isnan(cfmat)) ;
        cfmat = cfmat(:,indices) ; 
        cfmat = cfmat(indices,:) ; 
        vert = fieldnames_beh(indices) ; 
        horiz = vert' ; 
        cfmat_percent_fixed = cfmat ; 
        figure(figure_indicator)
        heatmap(horiz,vert,cfmat_percent_fixed,'XLabel','Actual Class','YLabel','Predicted Class','CellLabelFormat','%.1g') ; 
    
    case 2
        cfmat = cfmat_rawnum ; 
        indices = all(~isnan(cfmat));        
        cfmat = cfmat(:,indices) ; 
        cfmat = cfmat(indices,:) ; 
        vert = fieldnames_beh(indices) ; 
        horiz = vert' ; 
        cfmat_rawnum_fixed = cfmat ; 
        figure(figure_indicator)
        heatmap(horiz,vert,cfmat_rawnum_fixed,'ColorScaling','scaledcolumns', ...
            'XLabel','Actual Class','YLabel','Predicted Class') ; 

    case 3
        cfmat = cfmat_percent ; 
        indices = all(~isnan(cfmat)) ;
        cfmat = cfmat(:,indices) ; 
        indices2 = zeros(size(cfmat,1),1) ;
        for i = 1:length(indices2)
            if sum(cfmat(i,:)) == 0
                indices2(i) = 1; 
            end
        end
        indices3 = logical(~(indices2')) ; 
        cfmat = cfmat(indices3,:) ;
        vert = fieldnames_beh(indices3) ;
        horiz = fieldnames_beh(indices) ;
        cfmat_percent_fixed = cfmat ; 
        figure(figure_indicator)
        heatmap(horiz,vert,cfmat_percent_fixed,'XLabel','Actual Class','YLabel','Predicted Class','CellLabelFormat','%.1g')
        
    case 4
        cfmat = cfmat_rawnum ; 
        indices = all(~isnan(cfmat)) ;        
        indices2 = zeros(size(cfmat,1),1) ;
        for i = 1:length(indices2)
            if sum(cfmat(i,:)) == 0
                indices2(i) = 1; 
            end
        end
        indices3 = logical(~(indices2')) ; 
        vert = fieldnames_beh(indices3) ;
        horiz = fieldnames_beh(indices) ;        
        cfmat = cfmat(:,indices) ; 
        cfmat = cfmat(indices3,:) ; 
        cfmat_rawnum_fixed = cfmat ; 
        figure(figure_indicator)
        heatmap(horiz,vert,cfmat_rawnum_fixed,'ColorScaling','scaledcolumns', ...
            'XLabel','Actual Class','YLabel','Predicted Class') ;        
end

