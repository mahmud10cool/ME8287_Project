function [state, options, optchanged] = displayGeneration(options, state, flag)
    optchanged = false;

    if strcmp(flag, 'iter')
        fprintf('=== Generation %d ===\n', state.Generation);
    end
end