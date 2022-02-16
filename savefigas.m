%% Save figures from MATLAB into specified formats (such as colour eps).
%  f is the handle of the figure to save
%  filename is the name to give the file without the file extension
%  extension is the file format in which to save the figure
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

function savefigas(f,filename,extension)

% Check if the chosen extension requires a particular driver
switch extension
    case 'eps'
        driver = 'epsc';
    otherwise
        driver = '';
end

% Build the full filename, including extension
filename = strcat(filename,'.',extension);

% Call the saveas function with the appropriate number of arguments
if isempty(driver)
    saveas(f,filename);
else
    saveas(f,filename,driver);
end

end
