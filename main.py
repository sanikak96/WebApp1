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
   con = sql.connect("/home/sanika/proj/test.db")
   con.row_factory = sql.Row

   cur = con.cursor()
   cur.execute("select * from table3")

   rows = cur.fetchall();
   return render_template("datas.html",rows = rows)


@app.route('/predict/', methods=['GET', 'POST'])
def predict():



    if request.method == 'POST':

        seq = request.form['seq']
        with open("random.fasta", "w") as fp:
            fp.write(seq)

        pepdesc = PeptideDescriptor('/home/sanika/proj/random.fasta', 'eisenberg')  # use Eisenberg consensus scale
        globdesc = GlobalDescriptor('/home/sanika/proj/random.fasta')

        # --------------- Peptide Descriptor (AA scales) Calculations ---------------
        pepdesc.calculate_global()  # calculate global Eisenberg hydrophobicity
        pepdesc.calculate_moment(append=True)  # calculate Eisenberg hydrophobic moment

        # load other AA scales
        pepdesc.load_scale('gravy')  # load GRAVY scale
        pepdesc.calculate_global(append=True)  # calculate global GRAVY hydrophobicity
        pepdesc.calculate_moment(append=True)  # calculate GRAVY hydrophobic moment
        pepdesc.load_scale('z3')  # load old Z scale
        pepdesc.calculate_autocorr(1, append=True)  # calculate global Z scale (=window1 autocorrelation)

        # --------------- Global Descriptor Calculations ---------------
        globdesc.length()  # sequence length
        globdesc.boman_index(append=True)  # Boman index
        globdesc.aromaticity(append=True)  # global aromaticity
        globdesc.aliphatic_index(append=True)  # aliphatic index
        globdesc.instability_index(append=True)  # instability index
        globdesc.calculate_charge(ph=7.4, amide=False, append=True)  # net charge
        globdesc.calculate_MW(amide=False, append=True)  # molecular weight

        f1=pepdesc.descriptor
        f2=globdesc.descriptor
        result=np.concatenate((f2,f1),axis=1)

        clf=joblib.load('ml_model.pkl')
        pred=clf.predict(result)
	proba=clf.predict_proba(result).tocoo()
        mc=pred.tocoo()
        out=mc.col
        res=[]
	labels=['antiviral','antibacterial','antifungal']
        values=proba.data
        plt.pie(values,labels=labels,autopct='%.0f%%',shadow=True, radius=0.5)
        plt.savefig('/home/sanika/proj/pie_chart.jpg')
	
        figfile = BytesIO()
        plt.savefig(figfile, format='png')
        figfile.seek(0)
	figdata_png = base64.b64encode(figfile.getvalue()).decode('ascii')
	plt.close()
    
        for i in range(len(out)):
            if out[i]==0:
               res.append("antiviral")
            elif out[i]==1:
               res.append("antibacterial")
            else:
               res.append("antifungal")

	

        return render_template('seq.html', seq = res, result=figdata_png)

    return render_template('predictor.html')



@app.route('/predict/upload/', methods=['GET', 'POST'])

def upload():

    if request.method == 'POST':
        # This will be executed on POST request.
        upfile = request.files['file']
        if upfile and allowed_file(upfile.filename):
	    
            filename = secure_filename(upfile.filename)
            upfile.save(os.path.join(app.config['UPLOAD_FOLDER'], filename))
	    #return render_template('upload.html')
	    #flash("File uploaded", "success")
            #with open("/home/sanika/proj/uploads/aa.fasta") as f:
    		#lines = f.readlines()
    		#lines = [l for l in lines if "ROW" in l]


    	    #with open("/home/sanika/proj/uploads/out.fasta", "w") as f1:
        	#f1.writelines(lines)

	    #f = open(filename)
	    #prot_seq = ReadFasta(f)

	    with open(filename) as fasta_file:  # Will close handle cleanly
    		identifiers = []
    		sequence = []
    	   	for seq_record in SeqIO.parse(fasta_file, 'fasta'):  # (generator)
        		identifiers.append(seq_record.id)
        		sequence.append(seq_record.seq)
	    
	    pepdesc = PeptideDescriptor(filename, 'eisenberg')  # use Eisenberg consensus scale
            globdesc = GlobalDescriptor(filename)

        # --------------- Peptide Descriptor (AA scales) Calculations ---------------
            pepdesc.calculate_global()  # calculate global Eisenberg hydrophobicity
            pepdesc.calculate_moment(append=True)  # calculate Eisenberg hydrophobic moment

        # load other AA scales
            pepdesc.load_scale('gravy')  # load GRAVY scale
            pepdesc.calculate_global(append=True)  # calculate global GRAVY hydrophobicity
            pepdesc.calculate_moment(append=True)  # calculate GRAVY hydrophobic moment
            pepdesc.load_scale('z3')  # load old Z scale
            pepdesc.calculate_autocorr(1, append=True)  # calculate global Z scale (=window1 autocorrelation)

        # --------------- Global Descriptor Calculations ---------------
            globdesc.length()  # sequence length
            globdesc.boman_index(append=True)  # Boman index
            globdesc.aromaticity(append=True)  # global aromaticity
            globdesc.aliphatic_index(append=True)  # aliphatic index
            globdesc.instability_index(append=True)  # instability index
            globdesc.calculate_charge(ph=7.4, amide=False, append=True)  # net charge
            globdesc.calculate_MW(amide=False, append=True)  # molecular weight

            f1=pepdesc.descriptor
            f2=globdesc.descriptor
            result=np.concatenate((f2,f1),axis=1)
	    rs=[]
	    for i in range(len(result)):
	 	prt=np.reshape(result[i],(-1,14))
		clf=joblib.load('ml_model.pkl')
		pred=clf.predict(prt)
		out=pred.toarray()
	#print(clf.predict_proba(result))
		proba=clf.predict_proba(prt).tocoo()
		mc=pred.tocoo()
		out=mc.col
		res=[]
		for i in range(len(out)):
	    		if out[i]==0:
	       			res.append("antiviral")
	    		elif out[i]==1:
	       			res.append("antibacterial")
	    		else: 
	       			res.append("antifungal")
		rs.append(res)
	    a = []
	    for i in range(len(rs)):
      		a.append('-'.join(rs[i]))

	    df = pd.DataFrame(data={"id": identifiers, "sequence": sequence, "activity":a},columns=['id','sequence','activity'])
	    df.to_csv("result.csv", sep=',',index=False)

	   
	   
		
	   
	    


            

	  
            os.remove(os.path.join(app.config['UPLOAD_FOLDER'], filename))
	    
	    #return render_template('seq.html', seq = rs)
	    return render_template('up.html', mimetype="text/csv")
	    
		

    

            #flash("File uploaded: Thanks!", "success")
	else:
		error = "PLEASE CHECK THE FORMAT OF FILE TO UPLOAD"
		return render_template('upload.html', error=error)   

    # This will be executed on GET request.
    return render_template('predictor.html')



@app.route('/predict/download/')
def download():
    		return send_file('result.csv',
                     mimetype='text/csv',
                     attachment_filename='result.csv',
                     as_attachment=True)


@app.route('/help/')
def help():
    return render_template('help.html')


if __name__ == '__main__':
    app.run(host='127.0.0.1', port='1024', debug=True)


