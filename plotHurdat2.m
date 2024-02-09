function plotHurdat2(t_target, lon, lat, t, kt34_ne, kt34_se, kt34_sw, kt34_nw, kt50_ne, kt50_se, kt50_sw, kt50_nw)
%{
  === plotHurdat2.m ===
  Plots storm track, storm center, and 34kt and 50kt wind radii from the NHC HURDAT2 Best Track file
  at a specified time. This time need not match a time slice in the HURDAT2 file, as all relevant 
  quantities are interpolated. For use with readHurdat2.m.
  
  Inputs:
     t_target (datetime)   - time for which to generate plot; used to interpolate all other values
     lon (double)          - x coordinate of the storm
     lat (double)          - y coordinate of the storm
     t (datetime)          - times at which the storm was at locations (lon,lat) in best track file
     kt34_ne (double)      - northeast quadrant 34kt wind radius maximum extent (in nautical miles)
     kt34_se (double)      - southeast quadrant 34kt wind radius maximum extent (in nautical miles)
     kt34_sw (double)      - southwest quadrant 34kt wind radius maximum extent (in nautical miles)
     kt34_nw (double)      - northwest quadrant 34kt wind radius maximum extent (in nautical miles)
     kt50_ne (double)      - northeast quadrant 50kt wind radius maximum extent (in nautical miles)
     kt50_se (double)      - southeast quadrant 50kt wind radius maximum extent (in nautical miles)
     kt50_sw (double)      - southwest quadrant 50kt wind radius maximum extent (in nautical miles) 
%}
    
    % take note of initial hold state
    hold_on = ishold;
    % convert nautical miles into degrees of latitude
    kt34_ne_n = nm2deg(kt34_ne);
    kt34_se_s = nm2deg(kt34_se);
    kt34_sw_s = nm2deg(kt34_sw);
    kt34_nw_n = nm2deg(kt34_nw);
    kt50_ne_n = nm2deg(kt50_ne);
    kt50_se_s = nm2deg(kt50_se);
    kt50_sw_s = nm2deg(kt50_sw);
    kt50_nw_n = nm2deg(kt50_nw);
    % convert nautical miles into degrees of longitude
    num_times = size(t,1);
    kt34_ne_e = zeros(num_times,1); 
    kt34_se_e = zeros(num_times,1);
    kt34_sw_w = zeros(num_times,1);
    kt34_nw_w = zeros(num_times,1);
    kt50_ne_e = zeros(num_times,1);
    kt50_se_e = zeros(num_times,1);
    kt50_sw_w = zeros(num_times,1);
    kt50_nw_w = zeros(num_times,1);
    earth_radius = 3440.065; % in nautical miles
    earth_radius_hurricane_lat = earth_radius * cosd(lat);
    for j = 1:num_times
        kt34_ne_e(j) = nm2deg(kt34_ne(j),earth_radius_hurricane_lat(j));
        kt34_se_e(j) = nm2deg(kt34_se(j),earth_radius_hurricane_lat(j));
        kt34_sw_w(j) = nm2deg(kt34_sw(j),earth_radius_hurricane_lat(j));
        kt34_nw_w(j) = nm2deg(kt34_nw(j),earth_radius_hurricane_lat(j));
        kt50_ne_e(j) = nm2deg(kt50_ne(j),earth_radius_hurricane_lat(j));
        kt50_se_e(j) = nm2deg(kt50_se(j),earth_radius_hurricane_lat(j));
        kt50_sw_w(j) = nm2deg(kt50_sw(j),earth_radius_hurricane_lat(j));
        kt50_nw_w(j) = nm2deg(kt50_nw(j),earth_radius_hurricane_lat(j));
    end
    % plot track
    plot(lon,lat,'-k','LineWidth',1)
    hold on
    plot(lon,lat,'ko','LineWidth',1)
    % interpolate all track values at t_target
    lon_interp = interp1(t, lon, t_target);
    lat_interp = interp1(t, lat, t_target);
    kt34_ne_n_interp = interp1(t, kt34_ne_n, t_target);
    kt34_ne_e_interp = interp1(t, kt34_ne_e, t_target);
    kt34_se_s_interp = interp1(t, kt34_se_s, t_target);
    kt34_se_e_interp = interp1(t, kt34_se_e, t_target);
    kt34_sw_s_interp = interp1(t, kt34_sw_s, t_target);
    kt34_sw_w_interp = interp1(t, kt34_sw_w, t_target);
    kt34_nw_n_interp = interp1(t, kt34_nw_n, t_target);
    kt34_nw_w_interp = interp1(t, kt34_nw_w, t_target);
    kt50_ne_n_interp = interp1(t, kt50_ne_n, t_target);
    kt50_ne_e_interp = interp1(t, kt50_ne_e, t_target);
    kt50_se_s_interp = interp1(t, kt50_se_s, t_target);
    kt50_se_e_interp = interp1(t, kt50_se_e, t_target);
    kt50_sw_s_interp = interp1(t, kt50_sw_s, t_target);
    kt50_sw_w_interp = interp1(t, kt50_sw_w, t_target);
    kt50_nw_n_interp = interp1(t, kt50_nw_n, t_target);
    kt50_nw_w_interp = interp1(t, kt50_nw_w, t_target);
    % plot storm center
    scatter(lon_interp,lat_interp,'filled','MarkerFaceColor','k','MarkerEdgeColor','w')
    % plot northeast sections
    step = pi/50;
    theta = 0:step:.5*pi;
    x = kt34_ne_e_interp * cos(theta) + lon_interp;
    y = kt34_ne_n_interp * sin(theta) + lat_interp;
    plot(x,y,'Color','#C724B1','LineWidth',2);
    x = kt50_ne_e_interp * cos(theta) + lon_interp;
    y = kt50_ne_n_interp * sin(theta) + lat_interp;
    plot(x,y,'k','LineWidth',2);
    % plot southeast sections
    theta = 1.5*pi:step:2*pi;
    x = kt34_se_e_interp * cos(theta) + lon_interp;
    y = kt34_se_s_interp * sin(theta) + lat_interp;
    plot(x,y,'Color','#C724B1','LineWidth',2);
    x = kt50_se_e_interp * cos(theta) + lon_interp;
    y = kt50_se_s_interp * sin(theta) + lat_interp;
    plot(x,y,'k','LineWidth',2);
    % plot southwest sections
    theta = pi:step:1.5*pi;
    x = kt34_sw_w_interp * cos(theta) + lon_interp;
    y = kt34_sw_s_interp * sin(theta) + lat_interp;
    plot(x,y,'Color','#C724B1','LineWidth',2);
    x = kt50_sw_w_interp * cos(theta) + lon_interp;
    y = kt50_sw_s_interp * sin(theta) + lat_interp;
    plot(x,y,'k','LineWidth',2);
    % plot northwest sections
    theta = .5*pi:step:pi;
    x = kt34_nw_w_interp * cos(theta) + lon_interp;
    y = kt34_nw_n_interp * sin(theta) + lat_interp;
    plot(x,y,'Color','#C724B1','LineWidth',2);
    x = kt50_nw_w_interp * cos(theta) + lon_interp;
    y = kt50_nw_n_interp * sin(theta) + lat_interp;
    plot(x,y,'k','LineWidth',2);
    % plot line segments connecting quarter-circular regions
    plot([lon_interp + kt34_ne_e_interp lon_interp + kt34_se_e_interp],[lat_interp lat_interp],'Color','#C724B1','LineWidth',2)
    plot([lon_interp - kt34_sw_w_interp lon_interp - kt34_nw_w_interp],[lat_interp lat_interp],'Color','#C724B1','LineWidth',2)
    plot([lon_interp lon_interp],[lat_interp - kt34_se_s_interp lat_interp - kt34_sw_s_interp],'Color','#C724B1','LineWidth',2)
    plot([lon_interp lon_interp],[lat_interp + kt34_nw_n_interp lat_interp + kt34_ne_n_interp],'Color','#C724B1','LineWidth',2)
    plot([lon_interp + kt50_ne_e_interp lon_interp + kt50_se_e_interp],[lat_interp lat_interp],'k','LineWidth',2)
    plot([lon_interp - kt50_sw_w_interp lon_interp - kt50_nw_w_interp],[lat_interp lat_interp],'k','LineWidth',2)
    plot([lon_interp lon_interp],[lat_interp - kt50_se_s_interp lat_interp - kt50_sw_s_interp],'k','LineWidth',2)
    plot([lon_interp lon_interp],[lat_interp + kt50_nw_n_interp lat_interp + kt50_ne_n_interp],'k','LineWidth',2)
    % create legend
    h(1) = plot(NaN,NaN,'Color','#C724B1','LineWidth',2);
    h(2) = plot(NaN,NaN,'k','LineWidth',2);
    h(3) = plot(NaN,NaN,'-ok');
    h(4) = scatter(NaN,NaN,'filled','MarkerFaceColor','k','MarkerEdgeColor','w');
    legend(h,'NHC 34kt max extent','NHC 50kt max extent','NHC storm track','NHC storm center','Location','southeast','Color','w');
    % reset hold state if necessary
    if ~hold_on
        hold off
    end
end