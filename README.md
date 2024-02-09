All of this code also exists on the Matlab File Exchange. I'm just storing it on GitHub as well so I can use git to work on any future File Exchange submissions.  

Matlab File Exchange does integrate with GitHub, but the integration creates more problems than it solves. The choices are:
1. Matlab automatically downloads the newest code from the default branch, but it doesn't bump the version number or otherwise change your File Exchange page. This can result in multiple versions of the code having the same version number and such.
2. Matlab automatically updates your release information using the information from GitHub Releases. This is not good for me because the code I intend to store here is not part of one cohesive project, and I have no interest in making a bunch of tiny repos for each new file I want to contribute. If I try to store all of this disparate code here anyway, I'll end up having versioned updates for some files on the File Exchange with no actual changes. This will be irritating for users, and the problem will only compound itself with each new file I create.
3. Don't use the built-in integration and manually make changes to File Exchange after I make changes on GitHub.

I am going with option #3.
