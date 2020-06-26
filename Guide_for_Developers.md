# Guide for Developers

### Checkout branch, edit, test, and merge

- To checkout a new branch named 'Example_Name' based on 'master'

```$git checkout -b Example_Name master```

- To push local edits to remote for new checked out branch named 'Example_Name' for first time

```$git push --set-upstream origin Example_Name```

- To test code, merge feature branch 'Example_Name' to 'development'

```$git fetch origin```\
```$git checkout development```\
```$git merge origin/Example_Name ```

Then make test for the new edits on 'development'. 
If you're happy with your result, issue a pull require.
If not, checkout out your feacture branch, make changes as needed on feature branch, then merge to 'development' and test again.


