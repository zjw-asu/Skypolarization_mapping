function [zenith, azimuth, elevation] = computeSunPosition(observer, time)
%COMPUTESUNPOSITION  Compute the sun's position in the sky.
%
%   [zenith, azimuth, elevation] = computeSunPosition(observer, time)
%
%   Returns the sun's angular position for a given geographic location and
%   local clock time.  Azimuth is measured from geographic North, clockwise
%   (0 = North, 90 = East, 180 = South, 270 = West).
%   Zenith is 0 at the point directly overhead and 90 at the horizon.
%   Elevation = 90 - zenith; negative elevation means sun is below horizon.
%
%   Reference equations:
%     http://www.me.umn.edu/courses/me4131/LabManual/AppDSolarRadiation.pdf
%
%   INPUTS
%   ------
%   observer  struct with fields:
%     .latitude   degrees (N positive)
%     .longitude  degrees (E positive)
%     .gmtOffset  hours from UTC  (e.g. -7 = MST, +8 = CST)
%
%   time      struct with fields:
%     .year, .month, .day   calendar date
%     .hour    local clock time in decimal hours (e.g. 14.5 = 2:30 PM)
%     .dst     daylight saving offset in hours (0 or 1; default 0)

lat = observer.latitude;
lon = observer.longitude;
gmtOffset = observer.gmtOffset;

year  = time.year;
month = time.month;
day   = time.day;
hour  = time.hour;
dst   = 0;
if isfield(time, 'dst'), dst = time.dst; end

% ---- Day of year ---------------------------------------------------------
daysPerMonth = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];
isLeapYear   = (mod(year,4)==0 && mod(year,100)~=0) || mod(year,400)==0;
if isLeapYear, daysPerMonth(2) = 29; end

if month == 1
    dayOfYear = day;
else
    dayOfYear = day + sum(daysPerMonth(1 : month-1));
end
n = dayOfYear;

% ---- Solar declination (degrees) -----------------------------------------
% Angle between the sun's rays and the Earth's equatorial plane.
% Ranges from -23.45 deg (winter solstice) to +23.45 deg (summer solstice).
delta = 23.45 * sind(360/365 * (284 + n));

% ---- Equation of Time (hours) --------------------------------------------
% Correction for Earth's elliptical orbit and axial tilt.
B   = 360/365 * (n - 81);
EOT = 0.165*sind(2*B) - 0.126*cosd(B) - 0.025*sind(B);

% ---- Local Solar Time (hours) --------------------------------------------
% Converts local clock time to solar time using the observer's exact longitude.
standardMeridian = 15 * gmtOffset;                        % degrees
localSolarTime   = hour - (1/15)*(standardMeridian - lon) + EOT - dst;

% ---- Solar hour angle (degrees) ------------------------------------------
% 0 at solar noon; negative before noon (morning); positive after noon.
hourAngle = 15 * (localSolarTime - 12);

% ---- Sun elevation and zenith --------------------------------------------
sinElevation = sind(lat)*sind(delta) + cosd(lat)*cosd(delta)*cosd(hourAngle);
sinElevation = max(-1, min(1, sinElevation));  % clamp for numerical safety

elevation = asind(sinElevation);
zenith    = 90 - elevation;

% ---- Sun azimuth (from North, clockwise) ---------------------------------
% Derived from the sun's unit vector in the local East-North-Up frame:
%   East  component: -cos(delta)*sin(hourAngle)
%   North component:  sin(delta)*cos(lat) - cos(delta)*cos(hourAngle)*sin(lat)
% atan2(East, North) gives azimuth measured clockwise from North.
%
% Note: hourAngle is positive in the afternoon (sun moves westward),
%       so -sin(hourAngle) is positive in the morning (sun in the East).
sunEast  = -cosd(delta) .* sind(hourAngle);
sunNorth =  sind(delta) .* cosd(lat) - cosd(delta) .* cosd(hourAngle) .* sind(lat);

azimuth = mod(atan2d(sunEast, sunNorth), 360);

end
