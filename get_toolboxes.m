function get_toolboxes( varargin)
%GET_TOOLBOXES Download the File Exchange submissions this solver needs.
%
%   get_toolboxes           download anything missing
%   get_toolboxes( 'force') re-download even if already present
%
% The submissions are not redistributed with this repository - two of the six
% ship without any license text, so they are fetched at setup time instead.
% Each one is unpacked into both single_sheet/ and double_sheets/ under the
% folder name that the respective add_pathes.m expects.
%
% TriStream is patched on the way in: tsearch was removed from MATLAB, and the
% velocity arrays have to be cast to double. See the Troubleshooting section of
% the README.

force = any( strcmpi( varargin, 'force'));
root = fileparts( mfilename( 'fullpath'));

base_url = 'https://www.mathworks.com/matlabcentral/mlc-downloads/downloads/submissions';

% id, name, subfolder inside the zip, target under single_sheet, target under double_sheets
tb = { ...
  33731, 'Meshing a plate using four noded elements', 'Plate Mesh',      'Plate_Mesh', 'Meshing_a_plate_using_four_noded_elements'; ...
  23488, 'Sparse sub access',                         'sparsesubaccess', fullfile('FEM_sparse','FEM_sparse'), 'Sparse_sub_access'; ...
  47092, 'Vectorized multi-dimensional matrix mult.', 'mntimes',         'mntimes', 'Vectorized_Multi-Dimensional_Matrix_Multiplication'; ...
  11278, 'TriStream',                                 'TriStream',       'TriStream', 'TriStream'; ...
  22351, 'Quiver 5',                                  '',                'quiver5', 'quiver5'; ...
  15881, 'mmwrite',                                   '',                'mmwrite', 'mmwrite'; ...
};

dest_base = { fullfile( root, 'single_sheet', 'cores', 'ToolBoxes'), ...
              fullfile( root, 'double_sheets', 'ToolBoxes') };

tmp = tempname;
mkdir( tmp);
cleanup = onCleanup( @() rmdir( tmp, 's'));

fprintf( '\nFetching %d File Exchange submissions into ToolBoxes/\n\n', size( tb, 1));

for ii = 1:size( tb, 1)
    id = tb{ii,1}; label = tb{ii,2}; sub = tb{ii,3};
    targets = { fullfile( dest_base{1}, tb{ii,4}), fullfile( dest_base{2}, tb{ii,5}) };

    if ~force && all( cellfun( @isfolder, targets))
        fprintf( '  [skip] %-46s already present\n', label);
        continue
    end

    zip_path = fullfile( tmp, sprintf( '%d.zip', id));
    url = sprintf( '%s/%d/versions/1/download/zip', base_url, id);
    try
        websave( zip_path, url);
    catch err
        warning( 'get_toolboxes:download', ...
                 'Could not download %s (FX %d): %s\nFetch it by hand from https://www.mathworks.com/matlabcentral/fileexchange/%d', ...
                 label, id, err.message, id);
        continue
    end

    ex = fullfile( tmp, sprintf( '%d', id));
    unzip( zip_path, ex);
    src = ex;
    if ~isempty( sub)
        src = fullfile( ex, sub);
    end

    for jj = 1:numel( targets)
        if isfolder( targets{jj})
            rmdir( targets{jj}, 's');
        end
        mkdir( targets{jj});
        copyfile( fullfile( src, '*'), targets{jj});
        % keep the license alongside the code where the submission ships one
        lic = fullfile( ex, 'license.txt');
        if isfile( lic)
            copyfile( lic, targets{jj});
        end
    end

    if id == 11278
        patch_tristream( targets);
    end

    fprintf( '  [ ok ] %-46s FX %d\n', label, id);
end

fprintf( '\nDone. cd into single_sheet or double_sheets and run GUI (or exe).\n');
fprintf( 'mpgwrite is not fetched - it is only needed for movie_format = ''mpeg'',\n');
fprintf( 'and the default ''wmv'' path uses mmwrite.\n\n');

end


function patch_tristream( targets)
%PATCH_TRISTREAM tsearch was removed from MATLAB; cast the velocities to double.

old_cast = 'x=x(:)''; y=y(:)''; x0=x0(:)''; y0=y0(:)''; u=u(:)''; v=v(:)'';';
new_cast = 'x=x(:)''; y=y(:)''; x0=x0(:)''; y0=y0(:)''; u=double(u(:)''); v=double(v(:)'');';
old_srch = 'TRI = tsearch(x,y,tri'',Xbeg,Ybeg);';
new_srch = 'TRI = tsearchn([x.'' y.''], tri'',[Xbeg.'' Ybeg.'']);';

for jj = 1:numel( targets)
    f = fullfile( targets{jj}, 'TriStream.m');
    if ~isfile( f), continue, end
    s = fileread( f);
    s = strrep( s, old_cast, new_cast);
    s = strrep( s, old_srch, new_srch);
    fid = fopen( f, 'w');
    fwrite( fid, s);
    fclose( fid);
end

end
