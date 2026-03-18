class Person:

    def __init__(self, firstname, lastname):
        self.firstname = firstname
        self.lastname = lastname

    def printFullName(self):
        print(f"{self.firstname} {self.lastname}")

class Student(Person):

    def __init__(self, firstname, lastname, subject):
        super().__init__(firstname, lastname)
        self.subject = subject

    def printNameSubject(self):
        print(f"{self.firstname} {self.lastname}, {self.subject}")

class Teacher(Person):

    def __init__(self, firstname, lastname, course):
        super().__init__(firstname, lastname)
        self.course = course

    def printNameCourse(self):
        print(f"{self.firstname} {self.lastname}, {self.course}")