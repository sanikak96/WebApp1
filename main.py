from flask import Flask, render_template, request
import sqlite3 as sql


app = Flask(__name__)


@app.route('/')
def main():
    return render_template('index.html')


@app.route('/about/')
def about():
    return render_template('about.html')


@app.route('/dataset/')
def datas():
   con = sql.connect("/home/projectBE/mysite/test.db")
   con.row_factory = sql.Row

   cur = con.cursor()
   cur.execute("select * from table3")

   rows = cur.fetchall();
   return render_template("datas.html",rows = rows)


@app.route('/pred/', methods=['GET', 'POST'])
def pred():
    if request.method == "POST":
        seq = request.form['seq']
        return render_template('seq.html', seq=seq)

    return render_template('predictor.html')


@app.route('/help/')
def help():
    return render_template('help.html')


if __name__ == '__main__':
    app.run(host='127.0.0.1', port='1024')


