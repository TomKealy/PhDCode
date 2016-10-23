% A function which provides parameters description
function [description] = do_work2(frequency, mcs_mode)
    description = ['Frequency is ' num2str(frequency) 'MHz, mcs mode is ' num2str(mcs_mode)];
    result = frequency * mcs_mode;
    save(['Output/' num2str(frequency) '_' num2str(mcs_mode) '.mat'], 'result');
end