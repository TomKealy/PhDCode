% Simple function which does not provide parameters description
function [] = do_work1(value)
    save(['Output/' num2str(value) '.mat'], 'value');
end