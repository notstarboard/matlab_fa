function [t, record_id, status, lon, lat, kt_max, mbar_min, kt34_ne, kt34_se, kt34_sw, kt34_nw, kt50_ne, kt50_se,...
    kt50_sw, kt50_nw, kt64_ne, kt64_se, kt64_sw, kt64_nw, kt_max_rad] = readHurdat2(file, hurricane_id)
%{
  === readHurdat2.m ===
  Reads in all data for a specified storm from the NHC HURDAT2 Best Track file
  See https://www.nhc.noaa.gov/data/#hurdat for downloads and for information about the HURDAT2 format

  Inputs:
     file (string)         - name of your HURDAT2 file (+ filepath if it's not in the current folder)
     hurricane_id (string) - storm/basin ID of the storm you're interested in (e.g. 'AL142018' for Hurricane Michael 2018)
  Outputs: 
     t (datetime)          - time at which the storm was at location (lon,lat)
     record_id (string)    - record identifier (see HURDAT2 spec for details)
     status (string)       - status of system (see HURDAT2 spec for details)
     lon (double)          - x coordinate of the storm
     lat (double)          - y coordinate of the storm
     kt_max (double)       - maximum sustained wind (in knots)
     mbar_min (double)     - minimum pressure (in millibars)
     kt34_ne (double)      - northeast quadrant 34kt wind radius maximum extent (in nautical miles)
     kt34_se (double)      - southeast quadrant 34kt wind radius maximum extent (in nautical miles)
     kt34_sw (double)      - southwest quadrant 34kt wind radius maximum extent (in nautical miles)
     kt34_nw (double)      - northwest quadrant 34kt wind radius maximum extent (in nautical miles)
     kt50_ne (double)      - northeast quadrant 50kt wind radius maximum extent (in nautical miles)
     kt50_se (double)      - southeast quadrant 50kt wind radius maximum extent (in nautical miles)
     kt50_sw (double)      - southwest quadrant 50kt wind radius maximum extent (in nautical miles)
     kt50_nw (double)      - northwest quadrant 50kt wind radius maximum extent (in nautical miles)
     kt64_ne (double)      - northeast quadrant 64kt wind radius maximum extent (in nautical miles)
     kt64_se (double)      - southeast quadrant 64kt wind radius maximum extent (in nautical miles)
     kt64_sw (double)      - southwest quadrant 64kt wind radius maximum extent (in nautical miles)
     kt64_nw (double)      - northwest quadrant 64kt wind radius maximum extent (in nautical miles)
     kt_max_rad (double)   - radius of maximum wind (in nautical miles)
%}

    % read in file
    track = readtable(file,'Format','%s%s%s%s%s%s%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f');
    % determine relevant rows for the specified hurricane
    exact_match_mask = strcmp(track{:,1}, hurricane_id);
    lower_bound = find(exact_match_mask) + 1; % find should only return one row
    j = lower_bound;
    upper_bound = 0;
    while j < size(track,1) && upper_bound == 0
       if isnan(str2double(track{j,1}))
           upper_bound = j - 1;
       end
       j = j + 1;
    end
    if upper_bound == 0
        upper_bound = size(track,1);
    end
    % get time
    hours = string(track{lower_bound:upper_bound,2});
    hours(hours == "0") = "0000";
    t = datetime(strcat(string(track{lower_bound:upper_bound,1})," ",hours),'InputFormat','uuuuMMdd HHmm');
    % get record identifier and status
    record_id = string(track{lower_bound:upper_bound,3});
    status = string(track{lower_bound:upper_bound,4});
    % get coordinates
    lat_raw = string(track{lower_bound:upper_bound,5});
    lat_raw(strlength(lat_raw) == 4) = strcat("0", lat_raw(strlength(lat_raw) == 4));
    lon_raw = string(track{lower_bound:upper_bound,6});
    lon_raw(strlength(lon_raw) == 4) = strcat("00", lon_raw(strlength(lon_raw) == 4));
    lon_raw(strlength(lon_raw) == 5) = strcat("0", lon_raw(strlength(lon_raw) == 5));
    lat = zeros(length(lat_raw),1);
    lon = zeros(length(lon_raw),1);
    % convert N, S, E, W into +/-N and +/-E
    for j=1:length(lon_raw)
        if lat_raw(j) == "" 
            lat(j) = double("");
        elseif extract(lat_raw(j),5) == 'S' 
            lat(j) = -double(extractBetween(lat_raw(j),1,4));
        else
            lat(j) = double(extractBetween(lat_raw(j),1,4));        
        end
        if lon_raw(j) == ""
            lon(j) = double("");
        elseif extract(lon_raw(j),6) == 'W'
            lon(j) = -double(extractBetween(lon_raw(j),1,5));
        else
            lon(j) = double(extractBetween(lon_raw(j),1,5));
        end
    end
    % get max wind
    kt_max = track{lower_bound:upper_bound,7};
    % get min pressure
    mbar_min = track{lower_bound:upper_bound,8};
    % get 34kt, 50kt, and 64kt radii
    kt34_ne = track{lower_bound:upper_bound,9};
    kt34_se = track{lower_bound:upper_bound,10};
    kt34_sw = track{lower_bound:upper_bound,11};
    kt34_nw = track{lower_bound:upper_bound,12};
    kt50_ne = track{lower_bound:upper_bound,13};
    kt50_se = track{lower_bound:upper_bound,14};
    kt50_sw = track{lower_bound:upper_bound,15};
    kt50_nw = track{lower_bound:upper_bound,16};
    kt64_ne = track{lower_bound:upper_bound,17};
    kt64_se = track{lower_bound:upper_bound,18};
    kt64_sw = track{lower_bound:upper_bound,19};
    kt64_nw = track{lower_bound:upper_bound,20};
    % get radius of max wind
    kt_max_rad = track{lower_bound:upper_bound,21};
end