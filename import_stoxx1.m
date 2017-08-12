function [IBassets, IBliabilities, IBassetsBBG, assets, liabilities, equity] = import_stoxx1(year)
%enter desired year as string, from 2014 to 2016.
%% interbank stress tests import data

M=xlsread('STOXX_banks.xlsx',year,'C2:H45');
IBassets = M(:,5);
IBliabilities = M(:,6);
IBassetsBBG = M(:,4);
assets = M(:,3);
liabilities = M(:,2);
equity = M(:,1);

end