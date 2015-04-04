# vario
Community variograms and associated null models

# install

```
#install.packages(devtools)
library(devtools)
install_github("mcglinnlab/vario")
```
# How to contribute to this project

1) Fork the repo to your local GitHub account

2) Clone your forked version of the repo to your machine
```
git clone git@github.com:your_user_name/vario.git
```

3) Link your local repo back to the master branch on mcglinnlab
```
git remote add upstream git@github.com:mcglinnlab/vario.git
```

4) Create a new branch that will contain your changes using
```
git branch new_feature
```
where `new_feature` is whatever name you decide to call your branch - typically
a string that suggests what the branch is meant to accomplish

5) Checkout the branch using 
```
git checkout new_feature
```

6) Stage the files you have created to accomplish your task using 
```
git add mycode.R
```

7) Commit your changes using 
```
git commit -m "Add a function to test for autocorrelation"
```
where the information in quotes is your commit message

8) Push your changes using 
```
git push origin new_feature
``` 
this will put your local changes to your `new_feature` branch online at your
personal github page.

9) Submit a pull request Via your web-browser. Navigate to your forked copy of
`vario` and navigate to your branch `new_feature` once you're on that branch
there should be a button on the right called "Pull request" click it. Add any
relevant information. The other collaborators on the project will recive an email.

10) Allow for time for your collaborators to review your changes. If accepted 
your changes will be merged into the master branch on mcglinnlab. 

11) Now that your branch `new_feature` has been merged into master you should 
update your forked copy of the repo by fetching the update `upstream` branch using
```
git fetch upstream
```

12)  Merge those changes with your local copy
```
git merge upstream/master
```

13)  Push the merge back to your forked repo
```
git push origin master
```
Always repeat steps 11-13 before you begin work on the project in the future. This
ensures that your local copy does not become out-of-date with the main repository.

# Licence
Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
