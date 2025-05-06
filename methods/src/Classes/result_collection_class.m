classdef result_collection_class < handle
    properties
        result_cellarray = {};
        current_size = 0;
    end

    methods

        function obj = result_collection_class()
            
        end

        function obj = set_results(obj, result)
            obj.current_size = obj.current_size + 1;
            obj.result_cellarray{obj.current_size} = result;           
        end

        function result = return_results(obj, str_name)
            result = [];
            for ii = 1:obj.current_size
                if strcmp(str_name, obj.result_cellarray{ii}.name)
                    result =  obj.result_cellarray{ii};
                end
            end
        end

        function result = return_names(obj)
            result = {};
            for ii = 1:obj.current_size
                result{end+1} = obj.result_cellarray{ii}.name;
            end

        end

      
    end
end