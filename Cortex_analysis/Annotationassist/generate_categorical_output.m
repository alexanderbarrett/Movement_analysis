function [outputvector,outputfields,observed_behaviors] = generate_categorical_output(labeled_struct,totalframes)

outputvector = zeros(1,totalframes);
outputfields = fieldnames(labeled_struct);
observed_behaviors = cell(1,1);
for mm = 1:numel(outputfields)
   if strcmp(outputfields{mm},'Unannotated')
        outputvector(labeled_struct.(outputfields{mm})) = nan;
   else
    outputvector(labeled_struct.(outputfields{mm})) = mm;
   end
    
   if numel(outputvector(labeled_struct.(outputfields{mm})))
       observed_behaviors{1,size(observed_behaviors,2)+1} = (outputfields{mm});
   end
   
end

end