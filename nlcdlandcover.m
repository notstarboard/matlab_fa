function [land_rough, lon, lat] = nlcdlandcover(File, Wlon, Elon, Slat, Nlat, Export)
%{
    === nlcdlandcover.m ===
    This function retrieves land cover data from the National Land Cover Database's (NLCD) continental 
    US (CONUS) land cover files, i.e. nlcd_*_land_cover_l48_*.i*, and maps it to land roughess using
    FEMA Hazus land cover to z0 mappings. You must download the land cover files from the MRLC website: 
    https://www.mrlc.gov/data (link good as of 6/2022).

    ***IMPORTANT NOTE***: This function requires a lot of memory to run. You'll want 35 GiB RAM + swap at
    minimum. One could improve performance by using memmapfile(), but that approach would likely require
    more maintenance than relying on functions like readgeoraster(). It would also require specific 
    knowledge about how data is organized in the .ige file. Hopefully Matlab will add functionality to 
    readgeoraster() in the future to allow users to read subsets of rasters and thereby improve performance.
    
    Author: Josh Port
    Date: June 2022

    I took inspiration from the usgslandcover.m function by Bryan C. Lougheed. Thanks, Bryan!
    
    --Input parameters-----------------------------------------------------
    
    File is a string specifying the location of the NLCD land cover .img file.
    file, e.g. 'nlcd_2019_land_cover_l48_20210604.img'. The .ige file from the
    zip must also be present in the same folder.
    
    Wlon, Elon, Slat and Nlat are the decimal degree values of the
    western, eastern, southern and northern limits of the desired lat/lon
    window, whereby degrees N should be entered as positive, degrees S as
    negative, degrees W as negative and degrees E as positive. The nearest
    possible coordinates in the NLCD file will be selected.
    
    NOTE: Due to the map projection, cartesian rows and columns in the NLCD data are
          not parallel to lines of latitude and longitude. This function  
          does not interpolate the output onto a regular grid. It prioritizes data  
          accuracy over making the output easy to use.
    
    Export is a flag to determine whether to generate a NetCDF file
    containing lon, lat, and z0 for the specified region.
    
    --Output data----------------------------------------------------------
    
    land_rough: Matrix containing the gridded roughness length.
    
    lon, lat: Matrices of same dimensions as land_rough containing the
    longitude and latitude coordinates for each datapoint
    in land_rough. Coordinates are in degrees east and north respectively.
    
    --Example--------------------------------------------------------------
    
    Retrieve the roughness length data for Rhode Island:
    
    [land_rough, lon, lat] = ...
      nlcdlandcover('/home/user/Downloads/nlcd_2019_land_cover_l48_20210604.img', -72, -71, 41, 42.05, 1);
    
    --Roughness length mapping---------------------------------------------
    
    Land Cover to z0 mapping was retrieved from the Hazus Hurricane Model Technical Manual for
    Hazus 4.2 Service Pack 3 published in March 2021. See Tables 4-13 and 4-14.

    https://www.fema.gov/sites/default/files/documents/fema_hazus-hurricane-technical-manual-4.2.3_0.pdf 
    * Link good as of 6/2022.
    https://www.fema.gov/flood-maps/tools-resources/flood-map-products/hazus/user-technical-manuals
    * Try here if the first link fails. A new version could be out.

    ***IMPORTANT NOTE***: Roughness mappings may be less accurate at the shoreline, as any points
    outside of the state boundaries defined in usastatehi.shp will fall back to the Median land
    cover to roughness length mapping. If you have access to more precise state boundary data 
    and/or know of a clever way of padding sea-facing state borders such that all land points in
    the land cover file are contained within them but state borders still never overlap, this
    would be a worthy customization to make.
%}

    tic
    % load raster file
    disp("Reading raster file")
    [cover,metadata] = readgeoraster(File);
    disp("Cropping to region of interest")
    % get coordinates of corners of region of interest in meters
    [x_nw,y_nw] = projfwd(metadata.ProjectedCRS,Nlat,Wlon);
    [x_ne,y_ne] = projfwd(metadata.ProjectedCRS,Nlat,Elon);
    [x_se,y_se] = projfwd(metadata.ProjectedCRS,Slat,Elon);
    [x_sw,y_sw] = projfwd(metadata.ProjectedCRS,Slat,Wlon);
    % define raster region
    w_bound_m = metadata.XWorldLimits(1);
    e_bound_m = metadata.XWorldLimits(2);
    s_bound_m = metadata.YWorldLimits(1);
    n_bound_m = metadata.YWorldLimits(2);
    x_step = metadata.CellExtentInWorldX;
    y_step = metadata.CellExtentInWorldY;
    % generate full list of x and y values
    x_vals=w_bound_m+x_step/2:x_step:e_bound_m-x_step/2;
    y_vals=n_bound_m-y_step/2:-y_step:s_bound_m+y_step/2;
    % find col and row indices at nearest to each corner of desired window
    [~, nw_col] = min(abs(x_vals - (x_nw+1e-100)));
    [~, nw_row] = min(abs(y_vals - (y_nw+1e-100)));
    [~, ne_col] = min(abs(x_vals - (x_ne+1e-100)));
    [~, ne_row] = min(abs(y_vals - (y_ne+1e-100))); 
    [~, se_col] = min(abs(x_vals - (x_se-1e-100)));
    [~, se_row] = min(abs(y_vals - (y_se-1e-100)));
    [~, sw_col] = min(abs(x_vals - (x_sw-1e-100)));
    [~, sw_row] = min(abs(y_vals - (y_sw-1e-100)));
    % use these indices to set region's bounds
    start_row = min([nw_row ne_row]);
    start_col = min([nw_col sw_col]);
    end_row = max([sw_row se_row]);
    end_col = max([ne_col se_col]);
    % pull in land cover values for region of interest
    land_cover = double(cover(start_row:end_row,start_col:end_col)); % cast because uint8 does not support NaN
    clear cover
    % convert from cartesian to lat/lon and create grid
    [x_grid,y_grid] = meshgrid(x_vals(start_col:end_col),y_vals(start_row:end_row));
    clear x_vals y_vals
    [lat,lon] = projinv(metadata.ProjectedCRS,x_grid,y_grid);
    clear x_grid y_grid
    % set all points outside the region of interest to NaN
    mask = lon > Elon | lon < Wlon | lat > Nlat | lat < Slat;
    lon(mask) = NaN;
    lat(mask) = NaN;
    land_cover(mask) = NaN;
    % map output land cover array to z0
    %%% load relevant state borders
    coastal_states = {'Alabama';'Connecticut';'Delaware';'Florida';'Georgia';'Louisiana';'Maine';...
        'Maryland';'Massachusetts';'Mississippi';'New Hampshire';'New Jersey';'New York';...
        'North Carolina';'Pennsylvania';'Rhode Island';'South Carolina';'Texas';'Virginia'};
    states = readgeotable("usastatehi.shp");
    state_border = states(ismember(states.Name,coastal_states),1:2); % not perfect, but the best solution without requiring additional input files
    %%% determine which points fall within which mapped regions
    mapping_region = {'Alabama';'Connecticut';'Delaware';'Florida NE';'Florida Panhandle';'Florida SE';...
        'Florida W';'Georgia';'Louisiana';'Maine';'Maryland';'Massachusetts';'Mississippi';'New Hampshire';...
        'New Jersey';'New York Long Island';'New York Manhattan';'New York Rest';'North Carolina';'Pennsylvania';...
        'Rhode Island';'South Carolina';'Texas';'Virginia';'Median'};
    mapping_region_state = {'Alabama';'Connecticut';'Delaware';'Florida';'Florida';'Florida';...
        'Florida';'Georgia';'Louisiana';'Maine';'Maryland';'Massachusetts';'Mississippi';'New Hampshire';...
        'New Jersey';'New York';'New York';'New York';'North Carolina';'Pennsylvania';...
        'Rhode Island';'South Carolina';'Texas';'Virginia';'N/A'};
    region_to_state = containers.Map(mapping_region,mapping_region_state);
    %%%%% order mapping regions by approximate distance from the center of the input region to hugely improve performance in the mask loop
    mapping_region_approx_ctr_lat = [33.09;41.70;38.95;29.00;30.68;26.43;29.16;32.79;31.41;45.47;39.11;...
                                     42.39;32.86;43.78;40.30;40.86;40.79;42.90;35.75;40.86;41.72;34.01;...
                                     31.74;37.69;-32.00]; % Median location is across the world so it will always sort last  
    mapping_region_approx_ctr_lon = [-86.77;-72.58;-75.51;-81.20;-85.82;-80.80;-82.34;-83.26;-92.54;-69.12;-76.42;...
                                     -71.90;-89.67;-71.54;-74.45;-73.14;-73.96;-75.69;-79.38;-77.72;-71.47;-80.78;...
                                     -98.85;-78.38;115.87]; % LabelLat and LabelLon in usastatehi.shp could provide many of these, but not for regions within states, so kept it manual 
    num_regions = size(mapping_region,1);
    input_ctr_lon = zeros(num_regions,1) + (Wlon + Elon)/2;
    input_ctr_lat = zeros(num_regions,1) + (Slat + Nlat)/2;
    ctr_dist = distance(input_ctr_lat,input_ctr_lon,mapping_region_approx_ctr_lat,mapping_region_approx_ctr_lon);
    [~,sorted_order] = sort(ctr_dist,'ascend');
    %%%%% all of the following intrastate boundaries are approximate and manually defined; did not consult FEMA for their exact boundaries
    fl_ew_bound = -81.7;
    fl_ns_bound = 27.5;
    fl_panw_bound = -83.7;
    manhattan = geopolyshape([40.880 40.753 40.696 40.711 40.729 40.743 40.794 40.835 40.873 40.880],...
        [-73.930 -74.011 -74.027 -73.977 -73.971 -73.971 -73.929 -73.934 -73.910 -73.930]);
    long_island = geopolyshape([40.539 40.539 41.013 41.227 41.213 40.876 40.806 40.795 40.779 40.743 40.707 40.703 40.680 40.624 40.539],...
        [-74.049 -73.500 -71.805 -71.811 -72.349 -73.745 -73.784 -73.905 -73.936 -73.964 -73.977 -73.997 -74.022 -74.047 -74.049]);
    mask_no_data = land_cover == 0;
    mask_done = mask_no_data | isnan(land_cover);
    state_info = shaperead("usastatehi.shp");
    mask = cell(num_regions,1);
    for i = 1:num_regions
        ordered_idx = sorted_order(i);
        region_name = mapping_region(ordered_idx);
        state_name = values(region_to_state,region_name);
        disp(strcat("Generating logical mask for ",region_name))
        if region_name == "Median" % last possible loop
            mask{num_regions} = ~mask_done;
        elseif all(mask_done,'all')  % no more masks needed
            sorted_order = sorted_order([1:i end]);
            break
        else
            mask_region = false(size(lon));
            region_name = mapping_region(ordered_idx);
            state = state_border.Shape(state_border.Name == state_name);
            %%%%% filter down to area within state's bounding box to minimize points passed to isinterior()
            bounding_box = state_info(states.Name == state_name).BoundingBox;
            maybe_in_state = lon >= bounding_box(1,1) & lon <= bounding_box(2,1) & lat >= bounding_box(1,2) & lat <= bounding_box(2,2);
            if region_name == "New York Rest" % no in_region mask for NY Rest; Long Island & Manhattan have their own geopolyshapes and are much smaller
                points_to_check = ~mask_done & maybe_in_state & ~mask{sorted_order(i-2)} & ~mask{sorted_order(i-1)};
            elseif state_name == "Florida"
                if region_name == "Florida NE"
                    maybe_in_region = lon >= fl_ew_bound & lat >= fl_ns_bound;
                elseif region_name == "Florida Panhandle"
                    maybe_in_region = lon <= fl_panw_bound;
                elseif region_name == "Florida SE"
                    maybe_in_region = lon >= fl_ew_bound & lat < fl_ns_bound;
                elseif region_name == "Florida W"
                    maybe_in_region = lon > fl_panw_bound;
                end
                points_to_check = ~mask_done & maybe_in_state & maybe_in_region;
            else
                points_to_check = ~mask_done & maybe_in_state;
            end 
            %%%%% see if any points to check are within region
            if any(points_to_check,'all')
                reverseStr = '';
                for j = 1:size(lon,1) % necessary because isinterior() will run out of memory in R2022a if you pass it all points at once
                    to_map_lat = lat(j,:);
                    to_map_lon = lon(j,:);
                    gps_to_map = geopointshape(to_map_lat(points_to_check(j,:)),to_map_lon(points_to_check(j,:)));
                    if region_name == "New York Long Island"
                        mask_region(j,points_to_check(j,:)) = isinterior(long_island,gps_to_map);
                    elseif region_name == "New York Manhattan"
                        mask_region(j,points_to_check(j,:)) = isinterior(manhattan,gps_to_map);
                    else
                        mask_region(j,points_to_check(j,:)) = isinterior(state,gps_to_map);
                    end
                    % Display progress; thanks to Undocumented Matlab for this block of code
                    percentDone = 100 * j / size(lon,1);
                    msg = sprintf('Percent done: %3.1f', percentDone);
                    fprintf([reverseStr, msg]);
                    reverseStr = repmat(sprintf('\b'), 1, length(msg));
                end
                fprintf(reverseStr)
            end
            mask{ordered_idx} = mask_region;
            mask_done = mask_done | mask_region;
        end
    end
    clear mask_done maybe_in_state maybe_in_region
    %%% define state-specific land cover to z0 mappings
    land_cover_id = [11 12 21 22 23 24 31 32 41 42 43 52 71 81 82 90 95]; %51/72/73/74 are for Alaska only;
    lulc_to_z0 = zeros(size(mapping_region,1),size(land_cover_id,2));
    lulc_to_z0(1,:) = [0.003 0.012 0.09 0.33 0.53 0.35 0.09 0.09 0.9 0.9 0.9 0.12 0.04 0.05 0.05 0.9 0.09];
    lulc_to_z0(2,:) = [0.003 0.012 0.1 0.34 0.48 0.35 0.1 0.1 0.9 0.9 0.9 0.11 0.04 0.05 0.06 0.9 0.09];
    lulc_to_z0(3,:) = [0.003 0.012 0.1 0.35 0.53 0.38 0.1 0.1 0.9 0.9 0.9 0.12 0.04 0.05 0.05 0.9 0.1];
    lulc_to_z0(4,:) = [0.003 0.012 0.09 0.33 0.54 0.34 0.09 0.09 0.9 0.9 0.9 0.12 0.04 0.06 0.06 0.9 0.09];
    lulc_to_z0(5,:) = [0.003 0.012 0.09 0.33 0.53 0.35 0.09 0.09 0.9 0.9 0.9 0.12 0.04 0.05 0.05 0.9 0.09];    
    lulc_to_z0(6,:) = [0.003 0.012 0.09 0.34 0.57 0.38 0.09 0.09 0.9 0.9 0.9 0.12 0.04 0.05 0.06 0.9 0.04];
    lulc_to_z0(7,:) = [0.003 0.012 0.1 0.34 0.54 0.39 0.1 0.09 0.9 0.9 0.9 0.12 0.04 0.05 0.06 0.9 0.08];
    lulc_to_z0(8,:) = [0.003 0.012 0.1 0.35 0.57 0.35 0.1 0.1 0.9 0.9 0.9 0.11 0.05 0.05 0.06 0.9 0.1];
    lulc_to_z0(9,:) = [0.003 0.012 0.09 0.33 0.5 0.39 0.09 0.09 0.9 0.9 0.9 0.12 0.04 0.06 0.06 0.9 0.11];
    lulc_to_z0(10,:) = [0.003 0.012 0.12 0.29 0.53 0.31 0.12 0.12 0.9 0.9 0.9 0.11 0.04 0.06 0.06 0.9 0.09];
    lulc_to_z0(11,:) = [0.003 0.012 0.1 0.35 0.53 0.38 0.1 0.1 0.9 0.9 0.9 0.12 0.04 0.05 0.05 0.9 0.1];
    lulc_to_z0(12,:) = [0.003 0.012 0.12 0.36 0.59 0.51 0.12 0.12 0.9 0.9 0.9 0.13 0.04 0.05 0.06 0.9 0.09];
    lulc_to_z0(13,:) = [0.003 0.012 0.09 0.33 0.53 0.35 0.09 0.09 0.9 0.9 0.9 0.12 0.04 0.05 0.05 0.9 0.09];
    lulc_to_z0(14,:) = [0.003 0.012 0.12 0.29 0.53 0.31 0.12 0.12 0.9 0.9 0.9 0.11 0.04 0.06 0.06 0.9 0.09];
    lulc_to_z0(15,:) = [0.003 0.012 0.08 0.35 0.58 0.42 0.08 0.08 0.9 0.9 0.9 0.13 0.04 0.06 0.05 0.9 0.08];
    lulc_to_z0(16,:) = [0.003 0.012 0.09 0.42 0.62 0.44 0.09 0.09 0.9 0.9 0.9 0.13 0.05 0.05 0.06 0.9 0.1];
    lulc_to_z0(17,:) = [0.003 0.012 0.09 0.5 0.84 1.55 0.09 0.09 0.9 0.9 0.9 0.13 0.05 0.05 0.06 0.9 0.1];  
    lulc_to_z0(18,:) = [0.003 0.012 0.09 0.43 0.73 0.92 0.09 0.09 0.9 0.9 0.9 0.13 0.05 0.05 0.06 0.9 0.1];
    lulc_to_z0(19,:) = [0.003 0.012 0.1 0.35 0.55 0.33 0.1 0.1 0.9 0.9 0.9 0.1 0.04 0.06 0.06 0.9 0.1];
    lulc_to_z0(20,:) = [0.003 0.012 0.14 0.36 0.62 0.44 0.14 0.14 0.9 0.9 0.9 0.14 0.04 0.05 0.05 0.9 0.1];
    lulc_to_z0(21,:) = [0.003 0.012 0.1 0.34 0.48 0.35 0.1 0.1 0.9 0.9 0.9 0.11 0.05 0.05 0.06 0.9 0.09];
    lulc_to_z0(22,:) = [0.003 0.012 0.1 0.36 0.57 0.34 0.1 0.1 0.9 0.9 0.9 0.13 0.04 0.06 0.06 0.9 0.1];
    lulc_to_z0(23,:) = [0.003 0.012 0.1 0.35 0.55 0.44 0.1 0.1 0.9 0.9 0.9 0.1 0.04 0.04 0.06 0.9 0.1];
    lulc_to_z0(24,:) = [0.003 0.012 0.1 0.35 0.53 0.38 0.1 0.1 0.9 0.9 0.9 0.12 0.04 0.05 0.05 0.9 0.1];
    lulc_to_z0(25,:) = [0.003 0.012 0.1 0.35 0.54 0.38 0.1 0.1 0.9 0.9 0.9 0.12 0.04 0.05 0.06 0.9 0.1];
    masks_in_use = size(sorted_order,1);
    z0_map = cell(masks_in_use,1);
    for i = 1:masks_in_use
        z0_map{sorted_order(i)} = containers.Map(land_cover_id,lulc_to_z0(sorted_order(i),:)); 
    end
    %%% assign roughness values based on mask values 
    disp("Assigning land roughness values")
    land_rough = NaN(size(land_cover));
    land_rough(mask_no_data) = 0.003; % points with no data should only be over water
    clear mask_no_data
    for i = 1:masks_in_use
        land_rough(mask{sorted_order(i)}) = cell2mat(values(z0_map{sorted_order(i)},num2cell(land_cover(mask{sorted_order(i)}))));
    end
    % transpose from lat,lon to lon,lat
    disp("Transposing")
    lon = lon';
    lat = lat';
    land_rough = land_rough';
    % write results to NetCDF file
    if Export == 1
        disp("Writing results to NetCDF file")
        % define output file and remove any file with the same name if it exists
        nc_out = strcat("NLCD_z0_",num2str(Wlon),",",num2str(Nlat),"_",num2str(Elon),",",num2str(Slat),".nc");
        if exist(nc_out,'file')
            delete(nc_out)
        end
        % write roughness
        nccreate(nc_out,'land_rough',...
                'Dimensions',{'lon',size(land_rough,1),'lat',size(land_rough,2)},...
                'Format','classic');
        ncwrite(nc_out,'land_rough',land_rough)
        % write lon
        nccreate(nc_out,'lon',...
                'Dimensions',{'lon',size(land_rough,1),'lat',size(land_rough,2)},...
                'Format','classic');
        ncwrite(nc_out,'lon',lon)
        % write lat
        nccreate(nc_out,'lat',...
                'Dimensions',{'lon',size(land_rough,1),'lat',size(land_rough,2)},...
                'Format','classic');
        ncwrite(nc_out,'lat',lat)
        % write all attributes
        ncid=netcdf.open(nc_out,'WRITE');
        netcdf.reDef(ncid)
        varid=netcdf.inqVarID(ncid,'land_rough');
        netcdf.putAtt(ncid,varid,'long_name','Roughness Length')
        netcdf.putAtt(ncid,varid,'units','Meters')
        varid=netcdf.inqVarID(ncid,'lon');
        netcdf.putAtt(ncid,varid,'long_name','Longitude')
        netcdf.putAtt(ncid,varid,'units','Degrees East')
        varid=netcdf.inqVarID(ncid,'lat');
        netcdf.putAtt(ncid,varid,'long_name','Latitude')
        netcdf.putAtt(ncid,varid,'units','Degrees North')
        netcdf.close(ncid)
    end
    toc
end
