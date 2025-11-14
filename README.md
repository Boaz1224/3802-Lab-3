%Use this to test my NACA Generator code (Change 4406 to any 4-digit number)This works for NACAGen4 and NACAGen3

[xb, yb] = naca4_boundary('4406', 201, 1.0, 'cosine');

figure; plot(xb, yb, '-'); axis equal; grid on;
xlabel('x'); ylabel('y'); title('NACA 2412 boundary (clockwise, TE \rightarrow TE)');

%Sanity checks
fprintf('Total boundary points: %d\n', numel(xb));
fprintf('Start @ TE: (%.6f, %.6f)\n', xb(1), yb(1));
fprintf('End   @ TE: (%.6f, %.6f)\n', xb(end), yb(end));
