function testSpinConv(testtype, verbose)
%TESTSPINCONV   Comparing SpinConv with SpinCalc.
%
%    TESTSPINCONV(TESTTYPE) is equivalent to TESTSPINCONV(TESTTYPE, FALSE)
%
%    TESTSPINCONV(TESTTYPE, VERBOSE) compares either the output of SpinConv
%    and SpinCalc (if TESTTYPE == 'output') or their execution speed (if
%    TESTTYPE == 'speed').
%
%    76 different kinds of rotations (about x, y, z, and their
%    compositions) are converted from any kind of supported rotation
%    representation (DCM, EV, Q, EA121, EA232, EA313, EA131, EA212, EA323,
%    EA123, EA231, EA312, EA132, EA213, EA321) to any other kind. 
%
%    All the equations used by SpinConv and SpinCalc are tested, including
%    the four alternative kinds of equations used to convert from DCM to
%    quaternion. Both arrays with N=1 (single rotation) and arrays with N>1
%    (multiple rotations) are tested.
%
%    If VERBOSE == FALSE (default value), negligible differences smaller
%    than 1000 * EPS are not shown. If VERBOSE == TRUE, all differences are
%    shown. 
%
%    If TESTTYPE == 'speed', VERBOSE is ignored, N is increased up to 10000
%    by replicating and tiling the above mentioned set of rotations, and
%    only a selected subset of the above mentioned comparisons is
%    performed, to avoid redundancy.
%
%    Note: 
%        Output differences larger than 1000 * EPS are due to the fact
%        that, when SpinCalc performs a conversion to EV (Axis-angle
%        representation), in some cases it may return an angle MU larger
%        than 180 degrees, while SpinConv in that case equivalently
%        represents the same rotation as a rotation by an angle smaller
%        than 180 degrees (MU = 360 - MU) about an Euler vector pointing in
%        the opposite direction.
%        Output differences smaller than 1000 * EPS are due to the
%        different methods used by SpinConv and SpinCalc to convert angles
%        from degrees to radians and vice-versa. SpinConv uses the same
%        method used by MATLAB functions DEG2RAD and RAD2DEG (see
%        README.TXT). 
%
%    See also SpinConv, SpinCalc.

% 2013 February 25
% Paolo de Leva
% University "Foro Italico" 
% Rome, Italy
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Narg = nargin;
switch Narg
    case 0, testtype = 'output'; verbose=false;
    case 1, verbose=false;
    case 2
    otherwise, error( nargchk(0, 2, Narg) );
end
testtype = lower(testtype);
switch testtype
    case {'output', 'speed'}
        % Do nothing
    otherwise
        error(['Invalid value for first input argument ' ...
               '(use either ''ourput'' or ''speed'')'])
end

% Rotations about x, y, and z
R1 = rotations( pi/6);   % (3x3x3) by   30 degrees
R2 = rotations( pi*2/3); % (3x3x3) by  120 degrees
R3 = rotations(-pi/3);   % (3x3x3) by  -60 degrees
R4 = rotations(-pi*5/6); % (3x3x3) by -150 degrees
R = cat(3, R1, R2, R3, R4); % 3x3x12

% Permutation (x x x x   y y y y   z z z z)
neworder = [1 4 7 10  2 5 8 11  3 6 9 12];
R = R(:,:,neworder);
    
% Compositions of rotations about x, y, z
comb3 = combN([1 2 3 4], 3);
x = comb3(:,1);
y = comb3(:,2)+4;
z = comb3(:,3)+8;
for i = 4^3 : -1 : 1
    Rcomp(:,:,i) = R(:,:,z(i)) * R(:,:,y(i)) * R(:,:,x(i)); % 3x3x64
end
    
% Testing all kinds or rotations with all kinds of supported representations
clc
if strcmp(testtype, 'output')
    types  = {'DCM' 'EV' 'Q' ...
              'EA121' 'EA232' 'EA313' 'EA131' 'EA212' 'EA323' ...
              'EA123' 'EA231' 'EA312' 'EA132' 'EA213' 'EA321'};
    types2 = {'DCM' 'EV' 'Q' ...
              'EA123' 'EA231' 'EA312' 'EA132' 'EA213' 'EA321'};
    disp '1) Testing rotations about x'
    disp ' ', test_them([types2 {'EA313' 'EA212'}], R(:,:,1:4), testtype, verbose);
    disp ' '
    disp '2) Testing rotations about y'
    disp ' ', test_them([types2 {'EA121' 'EA323'}], R(:,:,5:8), testtype, verbose);
    disp ' '
    disp '3) Testing rotations about z'
    disp ' ', test_them([types2 {'EA232' 'EA131'}], R(:,:,9:12), testtype, verbose);
    disp ' '
    disp '4) Testing compositions (Rx*Ry*Rz)'
    disp ' ', test_them(types, Rcomp, testtype, verbose); % Compositions (Rx*Ry*Rz)
else
    clear R1 R2 R3 R4 R comb3 x y z
    idx = 1 : 64;
    idx = idx(ones(1,157), :);
    Rcomp = Rcomp(:,:,idx(:)');
    types  = {'DCM' 'EV' 'Q' 'EA121' 'EA123'};
    disp '1) Testing single rotation (N = 1)'
    disp ' '
    disp '   Conversion type         A             B           A-B           A/B'   
    test_it(types, Rcomp(:,:,1), testtype, verbose);
    disp ' ', disp '   A = SpinConv; B = SpinCalc; ms = milliseconds.'
    disp ' '
    disp '2) Testing multiple rotations (N = 10)'
    disp ' '
    disp '   Conversion type         A             B           A-B           A/B'   
    test_it(types, Rcomp(:,:,1:10), testtype, verbose);
    disp ' ', disp '   A = SpinConv; B = SpinCalc; ms = milliseconds.'
    disp ' '
    disp '3) Testing multiple rotations (N = 100)'
    disp ' '
    disp '   Conversion type         A             B           A-B           A/B'   
    test_it(types, Rcomp(:,:,1:100), testtype, verbose);
    disp ' ', disp '   A = SpinConv; B = SpinCalc; ms = milliseconds.'
    disp ' '
    disp '4) Testing multiple rotations (N = 1000)'
    disp ' '
    disp '   Conversion type         A             B           A-B           A/B'   
    test_it(types, Rcomp(:,:,1:1000), testtype, verbose);
    disp ' ', disp '   A = SpinConv; B = SpinCalc; ms = milliseconds.'
    disp ' '
    disp '5) Testing multiple rotations (N = 10000)'
    disp ' '
    disp '   Conversion type         A             B           A-B           A/B'   
    test_it(types, Rcomp(:,:,1:10000), testtype, verbose);
    disp ' ', disp '   A = SpinConv; B = SpinCalc; ms = milliseconds.'
end
    

function test_them(types, R, testtype, verbose)
% Test conversions from any supported representation type to any other.
    
    % STEP 1 - Single rotation (N=1)
    N = size(R, 3);
    verbose0 = false;
    for i = 1 : N
        test_it(types, R(:,:,i), testtype, verbose0);
    end

    % STEP 2 - Multiple rotations (N>1)
    test_it(types, R, testtype, verbose);


function test_it(types, R, testtype, verbose)
% Test conversions from any supported representation type to any other.

    % STEP 1 - Computing with SpinCalc all kinds of representations
    n = length(types);  
    rotation{1} = R;
    for j = n : -1 : 2
        conv = ['DCM' 'to' types{j}];
        rotation{j} = SpinCalc(conv, R, 10*eps, 1);
    end

    % STEP 2 - Comparing SpinConv with SpinCalc for all kinds of coversions
    for i = 1 : n
        for j = 1 : n
            if j~=i
                conv = [types{i} 'to' types{j}];
                switch testtype
                    case 'output', compare_output(conv, rotation{i}, verbose)
                    case 'speed',  compare_speed (conv, rotation{i})
                end
            end
        end
    end


function compare_speed(conv, rotation)

    % Comparing execution speed
    f0 = @() SpinCalc(conv, rotation, 10*eps, 0);
    f1 = @() SpinConv(conv, rotation, 10*eps, 0);
    t0 = timeit(f0);
    t1 = timeit(f1);
    ratio = 100 * (t1 / t0);
    delta = t1 - t0;
    blank = ' '; blanks = blank( ones(1,14-length(conv)) );
    conv = ['   ' conv blanks];
    disp ([conv sprintf('%11.3f ms',     t1*1000) ...
                sprintf('%11.3f ms',     t0*1000) ...
                sprintf('%+11.3f ms', delta*1000) ...
                sprintf('%10i', round(ratio)) '%'])


function compare_output(conv, rotation, verbose)

    % Comparing output
    out1 = SpinCalc(conv, rotation, 10*eps, 1);
    out2 = SpinConv(conv, rotation);
    diff = abs(out2-out1);
    maxdiff = max( diff(:) );
    if maxdiff > 1000*eps
        disp (['   ' conv ...
            ' - Diff. between SpinCalc and SpinConv: ' ...
            num2str(maxdiff) ' (see help)'])
    elseif verbose && maxdiff > eps
        disp (['   ' conv ...
            ' - Diff. between SpinCalc and SpinConv < 1000*EPS (see help)'])
    elseif verbose
        disp (['   ' conv ...
            ' - No diff. between SpinCalc and SpinConv'])
    end


function R = rotations(theta)
% Rotations of coordinate system relative to vector space, about x, y and z
% ("alias" transformations)
    s = sin(-theta);
    c = cos(-theta);
    Rx = [ 1  0  0;
           0  c -s;
           0  s  c ];
    Ry = [ c  0  s;
           0  1  0;
          -s  0  c ];
    Rz = [ c -s  0;
           s  c  0;
           0  0  1 ];
    R = cat(3, Rx, Ry, Rz);


function M = combN(V, N)
% N-compositions with repetition of elements taken from vector V
% See ALLCOMB, available on MATLAB Central File Exchange
% (http://www.mathworks.com/matlabcentral/fileexchange)

% List of all possible combinations of N elements
[Y{N:-1:1}] = ndgrid(V);
M = cat(N+1, Y{:});
% Reshape into matrix
M = reshape(M, [], N);