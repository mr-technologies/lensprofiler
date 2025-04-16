$ErrorActionPreference="Stop"
pushd "${PSScriptRoot}"
$failed=$true
try
{
	"Installing Python modules..."
	python -m pip install --upgrade pip
	pip install --upgrade Pillow

	"Application dependencies (except for OpenCV) has been successfully installed."
	$failed=$false
}
catch
{
	$_
}
finally
{
	popd
	if($failed)
	{
		Read-Host "Failed! Press Enter to exit"
	}
	else
	{
		Read-Host "Success! Press Enter to exit"
	}
}
