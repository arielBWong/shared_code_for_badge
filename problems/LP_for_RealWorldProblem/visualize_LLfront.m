function visualize_LLfront(XL, FL, FLC)

f1 = figure();

index =  any((FLC) > 1e-6, 2);
FL = FL(~index, :);
scatter(FL(:, 1), FL(:, 2), 30, 'red', "filled");
hold on;

pause(1);
close(f1);


end