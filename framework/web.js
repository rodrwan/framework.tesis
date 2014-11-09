var express    = require("express");
var app        = express();
var formidable = require('formidable');
var exec       = require('child_process').exec;
var _          = require('lodash');
var colors     = require('colors');

colors.setTheme({
  silly: 'rainbow',
  input: 'grey',
  verbose: 'cyan',
  prompt: 'grey',
  info: 'green',
  data: 'grey',
  help: 'cyan',
  warn: 'yellow',
  debug: 'blue',
  error: 'red'
});

var form = "<!DOCTYPE HTML><html><body>" +
"<form method='post' action='/upload' enctype='multipart/form-data'>" +
"<input type='file' name='image'/>" +
"<input type='submit' /></form>" +
"</body></html>";

app.set('port', process.env.PORT || 3000);
app.use(express.favicon());
app.use(express.logger('dev'));
app.use(express.methodOverride());
app.use(express.cookieParser('our super ultra mega secret key'));
app.use(express.bodyParser());
app.use(express.cookieSession()); // Express cookie session middleware
app.use(app.router);
// app.use(express.bodyParser({limit: '50mb'}))
// app.use(express.json({limit: '50mb'}));
// app.use(express.urlencoded({limit: '50mb'}));

// development only
if ('development' == app.get('env')) {
  app.use(express.errorHandler());
}

app.get('/', function (req, res){
  res.writeHead(200, {'Content-Type': 'text/html' });
  res.end(form);
});

var fs = require('fs');

app.post('/upload', function (req, res) {
  res.header("Access-Control-Allow-Origin", "*");
  res.header("Access-Control-Allow-Headers", "X-Requested-With");
  res.header('Access-Control-Allow-Methods', 'GET, POST, OPTIONS, PUT, PATCH, DELETE');
  res.header('Access-Control-Allow-Headers', 'X-Requested-With,content-type');
  res.header('Access-Control-Allow-Credentials', true);

  var globalPath = '/Users/rodrwan/Documents/Algorithms/tesis/voc-release5/framework/';
  console.log('Uploading file...'.info);
  var timestamp = new Date().getTime();
  var imageName = 'uploads/bottle-'+timestamp+'.jpg';
  try {
    var base64Data = req.body.image.replace(/^data:image\/png;base64,/, "");
    var threshold = parseFloat(req.body.threshold);
  } catch (error) {
    console.log(error.error);
    res.writeHead(200, {'Content-Type': 'text/json'});
    res.end(JSON.stringify(error));
  }
  console.log("Processing image: " + imageName);
  if (isNaN(threshold)) {
    var error = 'Invalid threshold!';
    res.writeHead(200, {'Content-Type': 'text/json'});
    res.end(JSON.stringify(error));
  } else {
    fs.writeFile(imageName, base64Data, 'base64', function (err) {
      console.log('File uploaded...'.info);
      var image = globalPath + imageName;
      var model = globalPath + '2011/bottle/bottle_final';

      console.log("Running detection in matlab with a threshold of: %d\n".info, threshold);
      var instruction = "matlab -nodesktop -nosplash -r \"app('{{image}}', '{{model}}', {{threshold}}); quit\"";
      var matlabInstruction = instruction.replace('{{image}}', image).replace('{{model}}', model).replace('{{threshold}}', threshold);

      var start = process.hrtime();

      function execute(command, callback) {
        exec(command, function (error, stdout, stderr) { callback(stdout); });
      };

      var elapsed_time = function(note) {
        var precision = 3; // 3 decimal places
        var elapsed = process.hrtime(start)[1] / 1000000; // divide by a million to get nano to milli
        console.log(process.hrtime(start)[0] + " s, " + elapsed.toFixed(precision) + " ms - " + note); // print message + time
        start = process.hrtime(); // reset the timer
      };

      execute(matlabInstruction, function (output) {
        console.log(output);
        var lines = _.filter(output.split("\n"), function (line) {
            return !line.match(/[a-zA-Z]+/) && line.match(/\d+/);
        });
        var points = [];
        _.each(lines, function (line) {
          line = line.replace(/\s+/g, ' ');
          columns = line.split(" ").slice(1, -2);
          points.push({
            x1: Math.ceil(columns[0]),
            y1: Math.ceil(columns[1]),
            x2: Math.ceil(columns[2]),
            y2: Math.ceil(columns[3])
          });
        });
        elapsed_time("Sending response!.");
        res.writeHead(200, {'Content-Type': 'text/json' });
        res.end(JSON.stringify(points));
      });
    });
  }
});

console.log("Sirviendo en localhost:8080".warn);
app.listen(8080)