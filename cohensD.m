%% Calculate's Cohen's d for effect sizes between two means
%  The calculations are in line with the reference:
%  J. Cohen, Statistical power analysis for the behavioural sciences, 2nd ed. 
%  Hillsdale, USA: Lawrence Erlbaum Associates, 1988. Calculations are 
%  currently implemented for two-tailed tests with independent samples only.
%  Copyright (C) 2022  Richard Hugh Moulton
%
%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with this program.  If not, see <https://www.gnu.org/licenses/>.

function d = cohensD(x1,x2)
    n1 = length(x1);                % size of sample 1
    n2 = length(x2);                % size of sample 2
    s1 = sum((x1-mean(x1)).^2);     % adjusted variance for sample 1
    s2 = sum((x2-mean(x2)).^2);     % adjusted variance for sample 2
    s  =  sqrt((s1+s2)/(n1+n2-2));  % pooled within sample estimate of the population standard deviation
    
    d = abs(mean(x1)-mean(x2))/s;   % Cohen's d
end
