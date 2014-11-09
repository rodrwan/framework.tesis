var exec = require('child_process').exec;
var _ = require('lodash');

function execute(command, callback){
        exec(command, function(error, stdout, stderr){ callback(stdout); });
};

var image = './uploads/bottle.jpg';
var model = '/Users/rodrwan/Documents/Algorithms/tesis/voc-release5/framework/2011/bottle/bottle_final';
var threshold = -0.5;
var instruction = "matlab -nodesktop -nosplash -r \"app('{{image}}', '{{model}}', {{threshold}}); quit\"";
var matlabInstruction = instruction.replace('{{image}}', image).replace('{{model}}', model).replace('{{threshold}}', threshold);

execute(matlabInstruction, function(output) {
    var lines = _.filter(output.split("\n"), function(line) {
        return !line.match(/[a-zA-Z]+/) && line.match(/\d+/);
    });
    var points = []
    _.each(lines, function(line) {
        line = line.replace(/\s+/g, ' ');
        columns = line.split(" ").slice(1, -2);
        points.push({x1: columns[0], y1: columns[1], x2: columns[2], y2: columns[3]})
    });
    console.log(points);
});

