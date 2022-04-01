import pandas as pd

df = pd.read_csv('del_under.tsv', sep='\t')
x = df[df.iloc[:, 0] == 'Start Date']
df.columns = ["Start Date","End Date","Response Type","IP Address","Progress","Duration (in seconds)","Finished","Recorded Date","Response ID","Recipient Last Name","Recipient First Name","Recipient Email","External Data Reference","Location Latitude","Location Longitude","Distribution Channel","User Language","Please select your TA","Please identify your gender:","Is English your native language?","What is your current Class Year","Which of the following best describes your reason for enrolling in this class:","To what extent do you agree with the following statements: - I regularly attend recitation","To what extent do you agree with the following statements: - My TA encourages me to actively participate by providing opportunities to: gather information, synthesize, analyze, solve problems or reflect.","To what extent do you agree with the following statements: - My TA checks the class understands a topic before moving on","To what extent do you agree with the following statements: - My TA communicates content in a clear manner","To what extent do you agree with the following statements: - My TA seems comfortable teaching over Zoom (select N/A if not applicable)","To what extent do you agree with the following statements: - My TA provides opportunities for informal feedback that do not involve grading? (feedback on non-graded items e.i. in class questions, discussions, example problems etc)","To what extent do you agree with the following statements: - There are opportunities for students to provide feedback for the TA and course","To what extent do you agree with the following statements: - I feel adequately supported in this course","To what extent do you agree with the following statements: - I know when my TA's office hours are","To what extent do you agree with the following statements: - My TA is available to provide help","To what extent do you agree with the following statements: - Feedback I receive from my TA is constructive","To what extent do you agree with the following statements: - Feedback I receive from my TA is friendly","To what extent do you agree with the following statements: - I am comfortable asking my TA for help","To what extent do you agree with the following statements: - My TA respects students of all race, gender, class, ability/disability, religion, language, geographic region, sexual orientation and political affiliation","To what extent do you agree with the following statements: - My TA makes me feel safe in my classroom (select N/A if you do not meet in person)","What are the things you appreciate the most about your TA? What aspects of their instruction/mentorship would you like us to celebrate?","Mention the areas of improvement, that you would like to highlight for your TA. Please be constructive, so that the instructor can understand the situation better and take appropriate actions.","For remote learning specifically, what does your TA do that is effective?","How could your TA improve their remote instruction?","Is there any feedback for other TA (TA other than your recitation that you have visited) of this SAME course ?"]
for n in set(df['Please select your TA']):
    # print(n)
    if not isinstance(n,str):
        continue
    sub = df[df['Please select your TA'] == n]
    # sub.columns = list(x)
    name = n.split(' - ')[-1]
    name = name.replace(' ', '_')
    sub = sub.iloc[:,22:]
    sub.to_csv('DelUndergradCourses/{}.csv'.format(name))
# DelUndergradCourses