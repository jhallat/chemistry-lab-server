from chemistry.mole import molar_mass
from flask import Flask, jsonify, request
from flask_cors import CORS

app = Flask(__name__)
CORS(app)

response = {'results' : [
    {
        'number': '1',
        'status': 'ERR',
        'result': 'not implemented'
    }
]}

#TODO eventually remove this endpoint
@app.route('/molar_mass/<string:formula>')
def get_molar_mass(formula):
    return str(molar_mass(formula))

@app.route('/run-commands', methods=['POST'])
def run_commands():
    request_data = request.get_json()
    commands = request_data['commands']
    results = []
    for command in commands:
        result = {}
        result['number'] = command['number']
        if command['command'] == 'molar_mass':
            result['status'] = 'OK'
            result['result'] = str(molar_mass(command['parameters']['formula']))
        else:
            result['status'] = 'ERR'
            result['result'] = command['command'] + ' not implemented'
        results.append(result)
    response = {'results': results}
    return jsonify(response)



if __name__ == '__main__':
    app.run(port=5000)
