%==========================================================================
% Items for SVO and JPE (=HF) task 
%==========================================================================

%% SVO items

% SL is slider version A
SL(:,:, 1) = [85 85 85 85 85 85 85 85 85
    85 76 68 59 50 41 33 24 15];

SL(:,:, 2) = [85 87 89 91 93 94 96 98 100
    15 19 24 28 33 37 41 46 50];

SL(:,:, 3) = [50  54 59 63 68 72 76 81 85
    100 98 96 94 93 91 89 87 85];

SL(:,:, 4) = [50  54 59 63 68 72 76 81 85
    100 89 79 68 58 47 36 26 15];

SL(:,:, 5) = [100 94 88 81 75 69 63 56 50
    50  56 63 69 75 81 88 94 100];

SL(:,:, 6) = [100 98 96 94 93 91 89 87 85
    50 54 59 63 68 72 76 81 85];

SL(:,:, 7) = [100 96 93 89 85 81 78 74 70
    50 56 63 69 75 81 88 94 100];

SL(:,:, 8) = [90 91 93 94 95 96 98 99 100
    100 99 98 96 95 94 93 91 90];

SL(:,:, 9) = [100 94 88 81 75 69 63 56 50
    70 74 78 81 85 89 93 96 100];

SL(:,:,10) = [100 99 98 96 95 94 93 91 90
    70 74 78 81 85 89 93 96 100];

SL(:,:,11) = [70 74 78 81 85 89 93 96 100
    100 96 93 89 85 81 78 74 70];

SL(:,:,12) = [50  56 63 69 75 81 88 94 100
    100 99 98 96 95 94 93 91 90];

SL(:,:,13) = [50  56 63 69 75 81 88 94 100
    100 94 88 81 75 69 63 56 50];

SL(:,:,14) = [100 96 93 89 85 81 78 74 70
    90 91 93 94 95 96 98 99 100];

SL(:,:,15) = [90  91 93 94 95 96 98 99 100
    100 94 88 81 75 69 63 56 50];

% SE is slider example
SE = [30 35 40 45 50 55 60 65 70
    80 70 60 50 40 30 20 10 0];

% SR is slider version B % 1st row = item number version A - 2nd row is
% corresponding item number version B 
SR_vec = [ 1, 6;
    2,  5;
    3,  4;
    4,  3;
    5,  2;
    6,  1;
    7, 15;
    8, 14;
    9, 13;
    10, 12;
    11, 11;
    12, 10;
    13,  9;
    14,  8;
    15,  7; ];

SR = SL( :, :, SR_vec(:,2) );



%% JPE Items 
% rating from 1 (least preferable) to 4 (most preferable)

% Haruno & Frith 36
% left is self / right is other
HF = [ 036 177
    023 164
    013 150
    006 134
    002 117
    000 100
    066 006
    083 002
    100 000
    117 002
    134 006
    150 013
    198 117
    194 134
    187 150
    177 164
    164 177
    150 187
    002 083
    006 066
    013 050
    023 036
    036 023
    050 013
    134 194
    117 198
    100 200
    083 198
    066 194
    050 187
    164 023
    177 036
    187 050
    194 066
    198 083
    200 100 ];