import matplotlib.pyplot as plt
import pandas as pd

x = ['a', 'b', 'c', 'd', 'e', 'f']
y = [2, 10, 2, 7, 6, 1]


def plot(x, y, title, figname):
    y, x = zip(*sorted(zip(y, x)))
    fig, ax = plt.subplots()
    plt.barh(x, y)
    plt.title(title)
    plt.tight_layout()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.savefig(figname)
    plt.show()


"""
What ways do you like to learn?
"""
x = ['Lecture', 'Interactive', 'By teaching:', 'Reading', 'Demos', 'Group work']
y = [2, 10, 2, 7, 6, 1]
plot(x, y, 'W1 Q1: What ways do you like to learn?', 'q1.png')

"""
What would my ideal classroom dynamics look like?
"""
x = ['Interactive', 'Discussion', 'Lecture', 'Group work', 'Informal', 'Friendly/Safe']
y = [13, 3, 1, 5, 2, 9]
plot(x, y, 'W2 Q1: What would my ideal classroom dynamics look like?', 'q2.png')

"""
What makes you want to participate?
"""
x = ['Interactive/collaborative\nclassroom environment', 'Judgment free', 'Questions to answer', 'My personal interest',
     'Openness/friendliness\nof instructor', 'Professor makes\na topic interesting', 'Positive feedback ']
y = [12, 4, 6, 8, 4, 7, 3]
plot(x, y, 'W3 Q1: What makes you want to participate?', 'q3.png')

"""
What new skills did I gain from this class, or what skills did I dig deeper into?
"""
x = ['How to create interactive sessions', 'Communicate effectively', 'Humor', 'Teaching tools', 'Lesson planning',
     'Classroom management', 'Importance of Openess', 'To explain things in many ways']
y = [5, 11, 1, 1, 4, 3, 1, 4]
plot(x, y, 'W4 Q2: What new skills did I gain from this\nclass, or what skills did I dig deeper into?', 'q4.png')

"""
Character Arc
"""
x = [1, 2, 3, 4, 5]
y = [.8, .4, .2, 1, .8]
fig, ax = plt.subplots()

plt.plot(x, y)
plt.yticks([0, 1], ['snarky', 'good'])
plt.xticks([1,2,3,4,5])
plt.xlabel('Week')
plt.ylabel('Quality of self reflection')
plt.title("The quality of Student X's self reflections")
plt.tight_layout()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.savefig('character_arc.png')
plt.show()

