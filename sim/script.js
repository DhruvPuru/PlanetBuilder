function parse_float(x)
{
	if (isNaN(parseFloat(x)) || !isFinite(x))
		throw "Invalid number: " + x;
	return +x;
}

function parse_integer(x)
{
	x = parse_float(x);
	if (Math.round(x) != x)
		throw "Invalid integer: " + x;
	return Math.round(x);
}

function process(data)
{
	var canvas = document.getElementById("canvas");
	var data = data.split("\n");
	for (var i = 0 ; i != data.length ; ++i)
		data[i] = data[i].split(",");
	if (data.length < 2)
		throw "Invalid data format: " + data.length;
	if (data[0].length != 4)
		throw "Invalid header format: " + data[0].length;
	var time    = parse_integer(data[0][0]);
	var energy  = data[0][1].trim();
	var mass    = data[0][2].trim();
	var refresh = parse_integer(data[0][3]);
	// find edge points for drawing
	var max_xy = 0.0;
	for (var i = 1 ; i != data.length ; ++i) {
		if (data[i].length != 9)
			throw "Invalid orbit format: " + data[i].length;
		var planet_r = parse_float(data[i][4]);
		var center_x = parse_float(data[i][5]);
		var center_y = parse_float(data[i][6]);
		var orbit_a  = parse_float(data[i][7]);
		var center_xy = Math.max(Math.abs(center_x), Math.abs(center_y));
		max_xy = Math.max(max_xy, center_xy + orbit_a + planet_r);
	}
	// clear the canvas compute drawing scale
	var ctx = canvas.getContext("2d");
	ctx.clearRect(0, 0, canvas.width, canvas.height);
	var middle_x = canvas.width  * 0.5;
	var middle_y = canvas.height * 0.5;
	var scale = (Math.min(middle_x, middle_y) - 15) / max_xy;
	// draw year and day
	ctx.font = "15px Arial";
	ctx.textAlign = "left";
	ctx.fillStyle = "black";
	ctx.fillText("Year: "   + Math.floor(1 + time / 365), 10, 15);
	ctx.fillText("Day: " +              (1 + time % 365), 10, 35);
	ctx.fillText("Energy: " +                     energy, 10, 55);
	ctx.fillText("Mass (max): " +                   mass, 10, 75);
	// draw planets & asteroids
	for (var i = 1 ; i != data.length ; ++i) {
		// parse parameters
		var name     =             data[i][0].trim();
		var color    =             data[i][1].trim();
		var planet_x = parse_float(data[i][2]);
		var planet_y = parse_float(data[i][3]);
		var planet_r = parse_float(data[i][4]);
		var center_x = parse_float(data[i][5]);
		var center_y = parse_float(data[i][6]);
		var orbit_a  = parse_float(data[i][7]);
		var orbit_b  = parse_float(data[i][8]);
		// rotation angle
		var orbit_minus_A = 0.0;
		if (orbit_a != orbit_b)
			orbit_minus_A = Math.PI - Math.atan2(center_y, center_x);
		// convert center coordinates
		center_x = middle_x + center_x * scale;
		center_y = middle_y - center_y * scale;
		planet_x = middle_x + planet_x * scale;
		planet_y = middle_y - planet_y * scale;
		// draw the name of the planet
		if (name) {
			ctx.font = "11px Arial";
			ctx.textAlign = "center";
			ctx.fillStyle = color;
			ctx.fillText(name, planet_x, planet_y + planet_r + 10.0);
		}
		// draw the ellipse of the orbit
		ctx.save();
		ctx.beginPath();
		ctx.translate(center_x, center_y);
		ctx.rotate(orbit_minus_A);
		ctx.scale(orbit_a * scale, orbit_b * scale);
		ctx.arc(0.0, 0.0, 1.0, 0.0, 2.0 * Math.PI);
		ctx.restore();
		ctx.lineWidth = 0.5;
		ctx.strokeStyle = color;
		ctx.stroke();
		// draw the planet
		ctx.beginPath();
		ctx.arc(planet_x, planet_y, planet_r, 0, 2.0 * Math.PI);
		ctx.fillStyle = color;
		ctx.fill();
	}
	return refresh;
}

var latest_version = -1;

function ajax(version, retries, timeout)
{
	var xhr = new XMLHttpRequest();
	xhr.onload = (function() {
		var refresh = -1;
		try {
			if (xhr.readyState != 4)
				throw "Incomplete HTTP request: " + xhr.readyState;
			if (xhr.status != 200)
				throw "Invalid HTTP status: " + xhr.status;
			refresh = process(xhr.responseText);
			if (latest_version < version)
				latest_version = version;
			else
				refresh = -1;
		} catch (message) { alert(message); }
		if (refresh >= 0)
			setTimeout(function() { ajax(version + 1, 10, 100); }, refresh);
	});
	xhr.onabort   = (function() { location.reload(true); });
	xhr.onerror   = (function() { location.reload(true); });
	xhr.ontimeout = (function() {
		if (version <= latest_version)
			console.log("AJAX timeout (version " + version + " <= " + latest_version + ")");
		else if (retries == 0)
			location.reload(true);
		else {
			console.log("AJAX timeout (version " + version + ", retries: " + retries + ")");
			ajax(version, retries - 1, timeout * 2);
		}
	});
	xhr.open("GET", "data.txt", true);
	xhr.responseType = "text";
	xhr.timeout = timeout;
	xhr.send();
}

ajax(0, 10, 100);
